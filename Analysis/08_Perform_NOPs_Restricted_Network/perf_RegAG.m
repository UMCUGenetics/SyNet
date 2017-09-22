function result = perf_RegAG(dataset_info, opt_info)

%% Load CV info
iTr = dataset_info.DatasetTr.UsedSample;
iTe = dataset_info.DatasetTe.UsedSample;

%% Loading Train GE data
fprintf('Loading train data: [%s]\n', dataset_info.DatasetTr.GEPath);
data = load(dataset_info.DatasetTr.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
xTr = data.Gene_Expression(iTr,:);
Gene_Name = data.Gene_Name;
lTr = (data.Patient_Label(iTr)==1)*2-1;
if ~isequal(lTr, dataset_info.DatasetTr.Patient_Label), error(); end

%% Loading Train GE data
fprintf('Loading test data: [%s]\n', dataset_info.DatasetTe.GEPath);
data = load(dataset_info.DatasetTe.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
xTe = data.Gene_Expression(iTe,:);
lTe = (data.Patient_Label(iTe)==1)*2-1;
if ~isequal(lTe, dataset_info.DatasetTe.Patient_Label), error(); end

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);
n_gene = size(zTr, 2);

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
result.Gene_Name = Gene_Name;
end


