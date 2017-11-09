function result = perf_TNMC(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
n_gene = size(xTr, 2);

%% Select top genes
fprintf('Evaluting [%d] individual genes.\n', n_gene);
pv_vec = ttest2Ex(xTr, lTr);

%% Selecting top genes
[SubNet_Score, scr_ind] = sort(-log10(pv_vec), 'Descend');
SubNet_List = num2cell(scr_ind)';

%% Traning the final model
if opt_info.K == 0
	fprintf('Training the adaptive model over [%d] features...\n', n_gene);
	[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
	[pred, opt_K] = nmc(xTr(:,scr_ind), lTr, xTe(:,scr_ind), struct('iCvPar', Fold_Index));
	fprintf('Best K is: %d\n', opt_K);
	result.opt_K = opt_K;
else
	n_feat = min([opt_info.K n_gene]);
	fprintf('Training the final model over [%d] features...\n', n_feat);
	xTr = xTr(:, scr_ind(1:n_feat));
	xTe = xTe(:, scr_ind(1:n_feat));
	pred = nmc(xTr, lTr, xTe);
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


