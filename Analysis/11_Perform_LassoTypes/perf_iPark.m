function result = perf_iPark(dataset_info, opt_info, Int_Type)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
Net_Adj = dataset_info.DatasetTr.Net_Adj;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
[~, n_gene] = size(xTr);
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);

%% Park Method
% hLst=unique(round(n_gene-logspace(log10(1),log10(n_gene), 10)+1));
height_lst = unique(round(linspace(1,n_gene, 10)));
% height_lst(height_lst==n_gene) = []; %% In this case, it is using all genes (i.e. classic method and not park)
n_height = numel(height_lst);

fprintf('Performing linkage algorithm.\n');
Net_Adj(1:n_gene+1:end) = 1; % Set diagonal to one
Linkage_Z = linkage(squareform(1-Net_Adj), 'average');
ClsInd_mat = zeros(n_height, n_gene);

ClsInd_cell = cell(n_height,1);
B_cell = cell(n_height,1);
fit_cell = cell(n_height,1);
AUC_lst = zeros(n_height,1);
for hi=1:n_height
    fprintf('[%02d/%02d]: Performing clustering, requesting for [%d] clusters, ', hi, n_height, height_lst(hi));
	cls_ind = cluster(Linkage_Z, 'maxclust', height_lst(hi));
    [~, ~, ClsInd_mat(hi,:)] = unique(cls_ind, 'Stable');
    n_snet = max(ClsInd_mat(hi,:));
	fprintf('got [%d] clusters. \n', n_snet);
	
	ClsInd_cell{hi} = cell(n_snet,1);
	for ci=1:n_snet
        ClsInd_cell{hi}{ci} = find(ClsInd_mat(hi,:)==ci);
	end
	[mTr, ~] = IntegGenes(ClsInd_cell{hi}, zTr, lTr, zTr, Fold_Index, Int_Type);
	
	%% Normalizing features
	mTr = zscore(mTr);
	
	%% Inner cross-validation
	%fprintf('Performing lasso on the training set.\n');
	result = LassoWithCV(@lassoEx, mTr, lTr, [], [], dataset_info.DatasetTr.iCvPar, opt_info.lasso_opt);
	B = result.B;
	fit = result.fit;
	clear result
	
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
n_snet = numel(opt_B);
SubNet_Trimmed = ClsInd_cell{opt_height};
[SubNet_Score, cls_sid] = sort(abs(opt_B), 'Descend');
SubNet_Trimmed = SubNet_Trimmed(cls_sid);

%% Generating meta-genes
fprintf('Generating Meta-features from %d sub-networks.\n', n_snet);
[mTr, mTe] = IntegGenes(SubNet_Trimmed, zTr, lTr, zTe, Fold_Index, Int_Type);

%% Normalizing features
n_feat = min([n_snet MAX_N_SUBNET]);
mTr = zscore(mTr(:, 1:n_feat));
mTe = zscore(mTe(:, 1:n_feat));
	
%% Testing the final performance
fprintf('Training the final model ... \n');
result = LassoWithCV(@lassoEx, mTr, lTr, [], [], dataset_info.DatasetTr.iCvPar, opt_info.lasso_opt);
opt_Bmat = result.B;
opt_fit = result.fit;
clear result

%% Evaluating the model
opt_B = opt_Bmat(:, opt_fit.IndexMinMSE);
g = mTr*opt_B;
tr_auc = getAUC(lTr, g, 50);

g = mTe*opt_B;
te_auc = getAUC(lTe, g, 50);
% [~, ~, ~, auc]=perfcurve(l(iTe), g, fit.posClass);

fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', te_auc*100);

result.B = opt_Bmat;
result.fit = opt_fit;
result.height = opt_height;
result.SubNet_Full = ClsInd_cell;
result.SubNet_List = SubNet_Trimmed;
result.SubNet_Score = SubNet_Score;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.Gene_Name = dataset_info.DatasetTr.Gene_Name;
end