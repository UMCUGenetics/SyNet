function result = perf_CFGLasso(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
Net_Adj = dataset_info.DatasetTr.Net_Adj;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
[n_Trsample, n_gene] = size(xTr);
[n_Tesample,      ~] = size(xTe);
if isfield(opt_info, 'MAX_SUBNET_SIZE')
	MAX_SUBNET_SIZE = opt_info.MAX_SUBNET_SIZE;
else
	MAX_SUBNET_SIZE = 10;
end
fprintf('[i] Using MAX SN Size = %d\n', MAX_SUBNET_SIZE);

%% Correcting the gene directions
[xTr, xTe] = CorrectGeneDirection(xTr, xTe, lTr);

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);

%% Generate Neighbor Sets
fprintf('Generating neighbor sets and subnetworks: \n');
Neig_cell = getNeighborsFromAdj(Net_Adj);

%% Generate Subnetworks
SubNet_Full = cell(n_gene, 1);
for gi=1:n_gene
	showprogress(gi, n_gene);
	SubNet_Full{gi} = getNetNeighborsBreadthFirst(Neig_cell, Neig_cell{gi}, MAX_SUBNET_SIZE, 1);
end

%% Remove empty subnetworks
SubNet_size = cellfun(@(x) numel(x), SubNet_Full);
SubNet_Trimmed = SubNet_Full(SubNet_size>0);
n_snet = numel(SubNet_Trimmed);
SubNet_Str = cell(n_snet, 1);
for si=1:n_snet
	SubNet_Str{si,1} = sprintf('%d,', sort(SubNet_Trimmed{si}));
end
[~, SubNet_UID] = unique(SubNet_Str, 'Stable');
SubNet_Trimmed = SubNet_Trimmed(SubNet_UID);
n_snet = numel(SubNet_Trimmed);
fprintf('Subnetworks are purified. [%d] are left.\n', n_snet);

%% Construct Composite features
n_cf = 4;
cTr = zeros(n_Trsample, n_snet*n_cf);
cTe = zeros(n_Tesample, n_snet*n_cf);
group_ind = ones(3, n_snet);
step = 1;
for si=1:n_snet
	group_ind(1, si) = step;
	for ci=1:numel(SubNet_Trimmed{si})
		cTr(:, step) = zTr(:, SubNet_Trimmed{si}(ci));
		cTe(:, step) = zTe(:, SubNet_Trimmed{si}(ci));
		step = step + 1;
	end
	for ci=1:n_cf
		switch ci
			case 1
				func = @(lst)  min(lst,[], 2);
			case 2
				func = @(lst)  max(lst,[], 2);
			case 3
				func = @(lst) mean(lst, 2);
			case 4
				func = @(lst)  std(lst,0, 2);
		end
		cTr(:, step) = func(zTr(:, SubNet_Trimmed{si}));
		cTe(:, step) = func(zTe(:, SubNet_Trimmed{si}));
		step = step + 1;
	end
	
	%% Int grp
	[cTr(:, step), cTe(:, step)] = IntegGenes(SubNet_Trimmed(si), zTr, lTr, zTe, Fold_Index, 'RI');
	step = step + 1;
	
	group_ind(2, si) = step - 1;
end
cTr(:, step:end) = [];
cTe(:, step:end) = [];

%% Normalization
fprintf('Normalizing data ...\n');
cTr = zscore(cTr);
cTe = zscore(cTe);

%% Train Group lasso
fprintf('Training group lasso over [%d] groups and [%d] features ...\n', n_snet, size(cTr,2));
lasso_opt = {'lassoType', 'sgt', 'CV', 5, 'iCvPar', Fold_Index, 'relTol', 5e-2, 'n_lC', 15, 'lC_ratio', 1e-2, ...
	'n_lG', 15, 'lG_ratio', 1e-2, 'group_ind', group_ind, 'lambdaType', 'no-grid', 'paroptions', statset('UseParallel',false), 'verbose', 1};
[opt_B, opt_fit] = lassoEx(cTr, lTr, lasso_opt{:});

%% Collect Subnet scores
SubNet_Score = zeros(n_snet, 1);
for si=1:n_snet
	SubNet_Score(si) = mean(abs(opt_B(group_ind(1:2,si), opt_fit.IndexMinMSE)));
end

%% Showing results
n_lam = size(opt_B, 2);
tr_auc_lam = zeros(1, n_lam);
te_auc_lam = zeros(1, n_lam);
for i=1:n_lam
	tr_auc_lam(i) = getAUC(lTr, cTr*opt_B(:,i), 50);
	te_auc_lam(i) = getAUC(lTe, cTe*opt_B(:,i), 50);
end
fprintf('    Index: '); fprintf('%5d  ', 1:n_lam); fprintf('\n');
fprintf('Train AUC: '); fprintf('%0.3f, ', 1-opt_fit.MSE); fprintf('\n');
fprintf(' Test AUC: '); fprintf('%0.3f, ', te_auc_lam); fprintf('\n');
fprintf('Optimal lam is: [%d]\n', opt_fit.IndexMinMSE);

%% Plot results
%{
plot(tr_auc_lam); hold on; plot(te_auc_lam); plot(1-opt_fit.MSE); plot(opt_fit.IndexMinMSE,0.7, '*')
LC_lst = unique(opt_fit.Lambda);
LG_lst = unique(opt_fit.LambdaG);
tr_img = [];
te_img = [];
for li=1:numel(opt_fit.Lambda)
	i = find(ismember(LC_lst, opt_fit.Lambda(li)));
	j = find(ismember(LG_lst, opt_fit.LambdaG(li)));
% 	tr_img(i,j) = 1 - opt_fit.MSE(li);
	te_img(i,j) = te_auc_lam(li);
end
colormap(jet);
subplot(1,2,1);
imagesc(tr_img);
subplot(1,2,2);
imagesc(te_img);
set(gca, 'XTick', 1:numel(LG_lst), 'XTickLabel', LG_lst, 'XTickLabelRotation', 30, 'YTick', 1:numel(LC_lst), 'YTickLabel', LC_lst);
xlabel('Group');
ylabel('Lambda');
%}

%% Storing results
result.B = opt_B;
result.fit = opt_fit;
result.tr_auc_lam = tr_auc_lam;
result.te_auc_lam = te_auc_lam;
result.tr_auc = tr_auc_lam(opt_fit.IndexMinMSE);
result.te_auc = te_auc_lam(opt_fit.IndexMinMSE);
result.SubNet_List = SubNet_Trimmed;
result.SubNet_Score = SubNet_Score;
result.SubNet_GInd = group_ind;
result.Gene_Name = dataset_info.DatasetTr.Gene_Name;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
end


