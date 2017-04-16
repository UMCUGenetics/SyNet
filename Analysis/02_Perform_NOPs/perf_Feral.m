function result = perf_Feral(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
nTr = dataset_info.DatasetTr.Net_Adj;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTe = dataset_info.DatasetTe.Patient_Label;
[n_TrSample, n_gene] = size(xTr);
n_TeSample = size(xTe, 1);
MAX_SUBNET_SIZE = 7;
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Correcting the gene directions
gene_dir = corr(xTr, lTr);
for gi=1:n_gene
	if gene_dir(gi)<0
		xTr(:, gi) = -xTr(:, gi);
		xTe(:, gi) = -xTe(:, gi);
	end
end

%% Generate Neighbor Sets
fprintf('Generating neighbor sets and subnetworks: ');
Neig_cell = getNeighborsFromAdj(nTr);

%% Generate Subnetworks
SubNet_Full = cell(n_gene, 1);
for gi=1:n_gene
	showprogress(gi, n_gene);
	SubNet_Full{gi} = getNetNeighborsBreadthFirst(Neig_cell, Neig_cell{gi}, 10, 1);
end

%% Remove empty subnetworks
sn_size = cellfun(@(x) numel(x), SubNet_Full);
SubNet_Trimmed = SubNet_Full(sn_size>0);
n_snet = numel(SubNet_Trimmed);
fprintf('Subnetworks are purified. [%d] are left.\n', n_snet);

%% Trimming the subnetworks
fprintf('Sorting the subnetworks based on prediction score: ');
snet_auc = zeros(n_snet, 1);
for si=1:n_snet
	showprogress(si, n_snet, 20);
	SubNet_Trimmed{si} = SubNet_Trimmed{si}(1:min([numel(SubNet_Trimmed{si}) MAX_SUBNET_SIZE]));
	mTr = mean(xTr(:, SubNet_Trimmed{si}), 2);
	snet_auc(si,1) = measureAUC(mTr, lTr, 20);
	%snet_auc(si,2) = measureAUC(mTe, lTe, 20);
end
[snet_scr, snet_sid] = sort(snet_auc(:,1), 'descend');
SubNet_Trimmed = SubNet_Trimmed(snet_sid);

%% Generating Meta-features from Top SubNetworks
fprintf('Generating Meta-features from %d sub-networks.\n', n_snet);
mTr = zeros(n_TrSample, MAX_N_SUBNET);
mTe = zeros(n_TeSample, MAX_N_SUBNET);
for si=1:min([MAX_N_SUBNET n_snet])
    mTr(:, si) = mean(xTr(:, SubNet_Trimmed{si}), 2);
	mTe(:, si) = mean(xTe(:, SubNet_Trimmed{si}), 2);
end

%% Normalization
zTr = zscore(mTr);
zTe = zscore(mTe);
n_meta = size(zTr, 2);

%% Traning the final model
fprintf('Training the final model over [%d] meta-features...\n', n_meta);
[opt_B, opt_fit] = lassoEx(zTr, lTr, opt_info.lasso_opt{:}, 'iCvPar', dataset_info.DatasetTr.iCvPar);
fprintf('Final training is done. [%d] non-zero features identified.\n', sum(abs(opt_B(:, opt_fit.IndexMinMSE))>0));

%% Evaluating the model
vec_B = opt_B(:, opt_fit.IndexMinMSE);
g = zTr*vec_B;
tr_auc = getAUC(lTr, g, 50);

g = zTe*vec_B;
te_auc = getAUC(lTe, g, 50);

fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', te_auc*100);

result.B = opt_B;
result.fit = opt_fit;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.tr_mf = mTr;
result.te_mf = mTe;
result.SubNet_List = SubNet_Trimmed;
result.SubNet_Score = snet_scr;

% for i=1:size(zTr,2) %###
% 	tr_auc(i,1) = getAUC(lTr, zTr(:,i), 50);
% 	te_auc(i,1) = getAUC(lTe, zTe(:,i), 50);
% end
% plot(tr_auc, 'b'); hold on
% plot(te_auc, 'r');
end


