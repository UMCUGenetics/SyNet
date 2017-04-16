function result = perf_iPark(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
nTr = dataset_info.DatasetTr.Net_Adj;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTe = dataset_info.DatasetTe.Patient_Label;
[n_TrSample, n_gene] = size(xTr);
n_TeSample = size(xTe, 1);
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Park Method
% hLst=unique(round(n_gene-logspace(log10(1),log10(n_gene), 10)+1));
height_lst = unique(round(linspace(1,n_gene, 10)));
height_lst(height_lst==n_gene) = []; %% In this case, it is using all genes (i.e. classic method and not park)
n_height = numel(height_lst);

fprintf('Performing linkage algorithm.\n');
nTr(1:n_gene+1:end) = 1; % Set diagonal to one
Z = linkage(squareform(1-nTr), 'average');
ClsInd_mat = zeros(n_height, n_gene);

ClsInd_cell = cell(n_height,1);
B_cell = cell(n_height,1);
fit_cell = cell(n_height,1);
z_cell = cell(n_height,1);
AUC_lst = zeros(n_height,1);
for hi=1:n_height
    fprintf('[%02d/%02d]: Performing clustering, requesting for [%d] clusters, ', hi, n_height, height_lst(hi));
    [~, ~, ClsInd_mat(hi,:)] = unique(cluster(Z, 'maxclust', height_lst(hi)));
    n_cluster = max(ClsInd_mat(hi,:));
	fprintf('got [%d] clusters. \n', n_cluster);
	
	mTr = zeros(n_TrSample, n_cluster);
	ClsInd_cell{hi} = cell(n_cluster,1);
	for ci=1:n_cluster
        ClsInd_cell{hi}{ci} = find(ClsInd_mat(hi,:)==ci);
		mTr(:, ci) = mean(xTr(:, ClsInd_cell{hi}{ci}), 2);
	end
	
	%% Normalizing features
	zTr = zscore(mTr);
	z_cell{hi} = zTr;
	
	%% Inner cross-validation
	%fprintf('Performing lasso on the training set.\n');
	[B, fit] = lassoEx(zTr, lTr, opt_info.lasso_opt{:}, 'iCvPar', dataset_info.DatasetTr.iCvPar);
	
	B_cell{hi} = B;
	fit_cell{hi} = fit;
	AUC_lst(hi) = 1 - fit.MSE(fit.IndexMinMSE);
	fprintf('+ Crossvalidated AUC for [%04d] clusters: %0.2f%%\n', height_lst(hi), AUC_lst(hi)*100);
end

[~, opt_height] = max(AUC_lst);
fprintf('Optimal height is [%d] with [%0.2f%%] AUC.\n', height_lst(opt_height), AUC_lst(opt_height)*100);

%% Sorting subnetworks by their score
opt_fit = fit_cell{opt_height};
opt_Bmat = B_cell{opt_height};
opt_B = opt_Bmat(:, opt_fit.IndexMinMSE);
n_cluster = numel(opt_B);
cls_ind = ClsInd_cell{opt_height};
[cls_scr, cls_sid] = sort(abs(opt_B), 'descend');
cls_ind = cls_ind(cls_sid);

mTr = zeros(n_TrSample, MAX_N_SUBNET);
mTe = zeros(n_TeSample, MAX_N_SUBNET);
for ci=1:min([MAX_N_SUBNET n_cluster])
	mTr(:, ci) = mean(xTr(:, cls_ind{ci}), 2);
	mTe(:, ci) = mean(xTe(:, cls_ind{ci}), 2);
end

%% Normalizing features
zTr = zscore(mTr);
% if ~isequal(zTr, z_cell{opt_height}), error('Meta features are not consistant.\n'); end
zTe = zscore(mTe);
	
%% Testing the final performance
fprintf('Training the final model ... \n');
[opt_Bmat, opt_fit] = lassoEx(zTr, lTr, opt_info.lasso_opt{:}, 'iCvPar', dataset_info.DatasetTr.iCvPar);
fprintf('Final training is done. [%d] non-zero features identified.\n', sum(abs(opt_Bmat(:, opt_fit.IndexMinMSE))>0));

%% Evaluating the model
opt_B = opt_Bmat(:, opt_fit.IndexMinMSE);
g = zTr*opt_B;
tr_auc = getAUC(lTr, g, 50);

g = zTe*opt_B;
te_auc = getAUC(lTe, g, 50);
% [~, ~, ~, auc]=perfcurve(l(iTe), g, fit.posClass);

fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', te_auc*100);

result.B = opt_Bmat;
result.fit = opt_fit;
result.height = opt_height;
result.SubNet_Full = ClsInd_cell;
result.SubNet_Trimmed = cls_ind;
result.SubNet_Score = cls_scr;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.tr_mf = mTr;
result.te_mf = mTe;
end