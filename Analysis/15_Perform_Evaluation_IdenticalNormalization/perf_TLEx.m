function result = perf_TLEx(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
n_gene = size(xTr, 2);
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Select top genes
fprintf('Evaluting [%d] individual genes.\n', n_gene);
pv_vec = ttest2Ex(xTr, lTr);

%% Selecting top genes
[SubNet_Score, scr_ind] = sort(-log10(pv_vec), 'Descend');
SubNet_List = num2cell(scr_ind);
n_feat = min([MAX_N_SUBNET n_gene]);
xTr = xTr(:, scr_ind(1:n_feat));
xTe = xTe(:, scr_ind(1:n_feat));

%% Traning the final model
fprintf('Training the final model over [%d] features...\n', n_gene);
result = LassoWithCV(@lassoEx, xTr, lTr, xTe, lTe, dataset_info.DatasetTr.iCvPar, opt_info.lasso_opt);
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);

result.SubNet_List = SubNet_List;
result.SubNet_Score = SubNet_Score;
result.Gene_Name = dataset_info.DatasetTr.Gene_Name;
end


