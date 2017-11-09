function result = perf_TAgNMC(dataset_info, opt_info)

%% Load CV info
iTr = dataset_info.DatasetTr.UsedSample;
iTe = dataset_info.DatasetTe.UsedSample;
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Loading Train GE data
fprintf('Loading train data: [%s]\n', dataset_info.DatasetTr.GEPath);
data = load(dataset_info.DatasetTr.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
xTr = data.Gene_Expression(iTr,:);
Gene_Name = data.Gene_Name;
if isfield(opt_info, 'UseRndGene')
	fprintf('Using [%d] random genes ...\n', MAX_N_SUBNET);
	Rnd_gind = randperm(size(xTr,2), MAX_N_SUBNET);
	xTr = xTr(:, Rnd_gind);
	Gene_Name = Gene_Name(Rnd_gind);
end
lTr = (data.Patient_Label(iTr)==1)*2-1;
if ~isequal(lTr, dataset_info.DatasetTr.Patient_Label), error(); end
n_gene = size(xTr, 2);

%% Loading Train GE data
fprintf('Loading test data: [%s]\n', dataset_info.DatasetTe.GEPath);
data = load(dataset_info.DatasetTe.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
xTe = data.Gene_Expression(iTe,:);
if isfield(opt_info, 'UseRndGene')
	xTe = xTe(:, Rnd_gind);
end
lTe = (data.Patient_Label(iTe)==1)*2-1;
if ~isequal(lTe, dataset_info.DatasetTe.Patient_Label), error(); end

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);

%% Select top genes
fprintf('Evaluting [%d] individual genes.\n', n_gene);
pv_vec = ttest2Ex(zTr, lTr);

%% Selecting top genes
[SubNet_Score, scr_ind] = sort(-log10(pv_vec), 'Descend');
SubNet_List = num2cell(scr_ind)';
n_feat = min([MAX_N_SUBNET n_gene]);
zTr = zTr(:, scr_ind(1:n_feat));
zTe = zTe(:, scr_ind(1:n_feat));

%% Traning the final model
fprintf('Training the final model over [%d] features...\n', n_feat);
pred = nmc(zTr, lTr, zTe);

%% Saving results
result.tr_auc = 1;
result.te_auc = getAUC(lTe, pred, 50);
if isfield(opt_info, 'UseRndGene')
	result.Rnd_gind = Rnd_gind;
end
result.SubNet_List = SubNet_List;
result.SubNet_Score = SubNet_Score;
result.Gene_Name = Gene_Name;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
end


