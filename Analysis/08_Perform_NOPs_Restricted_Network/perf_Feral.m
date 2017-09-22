function result = perf_Feral(dataset_info, opt_info, Int_Type)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
Net_Adj = dataset_info.DatasetTr.Net_Adj;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
n_gene = size(xTr, 2);
if isfield(opt_info, 'MAX_SUBNET_SIZE')
	MAX_SUBNET_SIZE = opt_info.MAX_SUBNET_SIZE;
else
	MAX_SUBNET_SIZE = 5;
end
fprintf('[i] Using MAX SN Size = %d\n', MAX_SUBNET_SIZE);
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

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

%% Generating Meta-features from Top SubNetworks
fprintf('Generating Meta-features from %d sub-networks.\n', n_snet);
[mTr, mTe] = IntegGenes(SubNet_Trimmed, zTr, lTr, zTe, Fold_Index, Int_Type);

%% Sorting sub-networks
fprintf('Sorting sub-networks\n');
cr_net = corr(mTr, lTr, 'Type', 'Spearman');
[SubNet_Score, snet_sid] = sort(abs(cr_net), 'Descend');
SubNet_Trimmed = SubNet_Trimmed(snet_sid);
mTr = mTr(:, snet_sid);
mTe = mTe(:, snet_sid);

%% Normalization
n_feat = min([n_snet MAX_N_SUBNET]);
zTr = zscore(mTr(:,1:n_feat));
zTe = zscore(mTe(:,1:n_feat));
n_meta = size(zTr, 2);

%% Traning the final model
fprintf('Training the final model over [%d] meta-features...\n', n_meta);
result = LassoWithCV(@lassoEx, zTr, lTr, zTe, lTe, dataset_info.DatasetTr.iCvPar, opt_info.lasso_opt);
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);

%% Storing results
result.SubNet_List = SubNet_Trimmed;
result.SubNet_Score = SubNet_Score;
result.Gene_Name = dataset_info.DatasetTr.Gene_Name;
end


