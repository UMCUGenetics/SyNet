function result = perf_Regress(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
n_gene = size(xTr, 2);

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);

%% Traning the final model
fprintf('Training regression over [%d] features...\n', n_gene);
warning off
result.B = regress(lTr, zTr);
warning on
result.fit.IndexMinMSE = 1;

%% Evaluation
result.tr_auc = getAUC(lTr, zTr*result.B, 50);
result.te_auc = getAUC(lTe, zTe*result.B, 50);
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);

%% Saving the output
result.SubNet_List = num2cell(1:n_gene)';
result.SubNet_Score = (1:n_gene)';
result.Gene_Name = dataset_info.DatasetTr.Gene_Name;
end


