function result = perf_TTNMC(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
n_gene = size(xTr, 2);
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);

%% Select top genes
fprintf('Evaluting [%d] individual genes.\n', n_gene);
pv_mat = zeros(n_gene, 1);
for gi=1:n_gene
	[~, pv_mat(gi)] = ttest(zTr(:,gi), lTr);
end

%% Selecting top genes
[SubNet_Score, scr_ind] = sort(-log10(pv_mat), 'Descend');
SubNet_List = num2cell(scr_ind)';

%% Traning the final model
if isfield(opt_info, 'FindK')
	fprintf('Training the adaptive model over [%d] features...\n', n_gene);
	[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
	[pred, opt_K] = nmc(zTr(:,scr_ind), lTr, zTe(:,scr_ind), struct('iCvPar', Fold_Index, 'MAX_N_Feat', MAX_N_SUBNET));
	fprintf('Best K is: %d\n', opt_K);
	result.opt_K = opt_K;
else
	n_feat = min([MAX_N_SUBNET n_gene]);
	fprintf('Training the final model over [%d] features...\n', n_feat);
	zTr = zTr(:, scr_ind(1:n_feat));
	zTe = zTe(:, scr_ind(1:n_feat));
	pred = nmc(zTr, lTr, zTe);
	result.opt_K = n_feat;
end

%% Evaluating the model
tr_auc = 1;
te_auc = getAUC(lTe, pred, 50);

%% Saving results
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.SubNet_List = SubNet_List;
result.SubNet_Score = SubNet_Score;
result.Gene_Name = dataset_info.DatasetTr.Gene_Name;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
end


