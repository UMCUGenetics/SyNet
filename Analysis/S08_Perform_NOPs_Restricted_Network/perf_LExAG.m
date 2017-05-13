function result = perf_LExAG(dataset_info, opt_info)

%% Load CV info
iTr = dataset_info.DatasetTr.UsedSample;
iTe = dataset_info.DatasetTe.UsedSample;

%% Loading Train GE data
fprintf('Loading train data: [%s]\n', dataset_info.DatasetTr.GEPath);
data = load(dataset_info.DatasetTr.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
xTr = data.Gene_Expression(iTr,:);
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

%% Traning the final model
fprintf('Training the final model ...\n');
result = LassoWithCV(@lassoEx, zTr, lTr, zTe, lTe, dataset_info.DatasetTr.iCvPar, opt_info.lasso_opt);
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
result.Gene_Name = data.Gene_Name;
end


