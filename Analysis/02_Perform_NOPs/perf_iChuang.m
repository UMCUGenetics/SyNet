function result = perf_iChuang(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
nTr = dataset_info.DatasetTr.Net_Adj;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTe = dataset_info.DatasetTe.Patient_Label;
[n_TrSample, n_gene] = size(xTr);
n_TeSample = size(xTe, 1);
MAX_SUBNET_SIZE = 100;
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

GAMMA = 0.05; % According to Author's paper
n_bin = floor(log2(n_TrSample)+1); % According to Author's paper

%% Generate Neighbor Sets
fprintf('Generating neighbor sets and subnetworks: ');
Neig_cell = getNeighborsFromAdj(nTr);

%% Generate Subnetworks
SubNet_Full = cell(n_gene, 1);
SubNet_Trimmed = cell(n_gene,1);
SubNet_MI = zeros(n_gene, 1);
for gi=1:n_gene
	showprogress(gi, n_gene);
	SubNet_Full{gi} = getNetNeighborsBreadthFirst(Neig_cell, Neig_cell{gi}, MAX_SUBNET_SIZE, 1);
	[SubNet_Trimmed{gi}, SubNet_MI(gi)] = trimSubNetwork(SubNet_Full{gi}, xTr, lTr, GAMMA, n_bin);
end

%% Filter 1: Permute Expression values, P-value > 0.05
valid_subnet = filterSubNetworks_PermExprRndNode(SubNet_MI, 0.05, xTr, lTr, SubNet_Full, GAMMA, n_bin);
SubNet_Full = SubNet_Full(valid_subnet);
SubNet_Trimmed = SubNet_Trimmed(valid_subnet);
SubNet_MI = SubNet_MI(valid_subnet);

%% Filter 2: Start from same node, P-value > 0.05
valid_subnet = filterSubNetworks_PermExprSameNode(SubNet_MI, 0.05, xTr, lTr, SubNet_Full, GAMMA, n_bin);
SubNet_Full = SubNet_Full(valid_subnet);
SubNet_Trimmed = SubNet_Trimmed(valid_subnet);
SubNet_MI = SubNet_MI(valid_subnet);

%% Filter 3: Permuted labels, P-value > 0.00005
[valid_subnet, SubNet_PVal] = filterSubNetworks_PermLabels(SubNet_Trimmed, SubNet_MI, 0.00005, xTr, lTr, n_bin);
SubNet_Full = SubNet_Full(valid_subnet);
SubNet_Trimmed = SubNet_Trimmed(valid_subnet);
SubNet_MI = SubNet_MI(valid_subnet);
SubNet_PVal = SubNet_PVal(valid_subnet);

%% Retaining the unique subnetworks
fprintf('Selecting a unique set of sub-networks: ');
SubNet_Map = containers.Map();
n_snet = numel(SubNet_Trimmed);
valid_subnet = true(n_snet, 1);
for si=1:n_snet
    SubNet_Key = num2str(sort(SubNet_Trimmed{si}));
    if ~SubNet_Map.isKey(SubNet_Key)
        SubNet_Map(SubNet_Key) = si;
	else
		valid_subnet(si) = 0;
    end
end
SubNet_Full = SubNet_Full(valid_subnet);
SubNet_Trimmed = SubNet_Trimmed(valid_subnet);
SubNet_MI = SubNet_MI(valid_subnet);
SubNet_PVal = SubNet_PVal(valid_subnet);
fprintf('Filtering finished [%d --> %d].\n', n_snet, sum(valid_subnet));

%% Sorting SubNetworks
[~, sind] = sort(SubNet_PVal);
SubNet_Full = SubNet_Full(sind);
SubNet_Trimmed = SubNet_Trimmed(sind);
SubNet_MI = SubNet_MI(sind);
SubNet_PVal = SubNet_PVal(sind);

%% Meta feature generation
n_snet = numel(SubNet_Trimmed);
mTr = zeros(n_TrSample, MAX_N_SUBNET);
mTe = zeros(n_TeSample, MAX_N_SUBNET);
for fi=1:min([MAX_N_SUBNET n_snet])
    mTr(:,fi) = mean(xTr(:, SubNet_Trimmed{fi}),2);
	mTe(:,fi) = mean(xTe(:, SubNet_Trimmed{fi}),2);
end

%% Normalization
zTr = zscore(mTr);
zTe = zscore(mTe);
n_meta = size(zTr, 2);

%% Traning the final model
fprintf('Training the final model over [%d] meta-features...\n', n_meta);
old_AUC = 0;
for fi=1:n_meta
	[B, fit] = lassoEx(zTr(:, 1:fi), lTr, opt_info.lasso_opt{:}, 'iCvPar', dataset_info.DatasetTr.iCvPar);
	new_AUC = 1 - fit.MSE(fit.IndexMinMSE);
	fprintf('[%03d/%03d] Current AUC is: %0.2f%%\n', fi, n_meta, new_AUC*100);
	if old_AUC>new_AUC
		break;
	else
		opt_NF = fi;
		opt_B = B;
		opt_fit = fit;
	end
	old_AUC = new_AUC;
end
if isequal(fi, n_meta) && old_AUC<=new_AUC % If no break has happend
	opt_NF = fi;
	opt_B = B;
	opt_fit = fit;
end
fprintf('Final training is done. [%d] non-zero features identified.\n', sum(abs(opt_B(:, opt_fit.IndexMinMSE))>0));

%% Evaluating the model
vec_B = opt_B(:, opt_fit.IndexMinMSE);
g = zTr(:, 1:opt_NF)*vec_B;
tr_auc = getAUC(lTr, g, 50);

g = zTe(:, 1:opt_NF)*vec_B;
te_auc = getAUC(lTe, g, 50);

fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', te_auc*100);

result.B = opt_B;
result.fit = opt_fit;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.tr_mf = mTr;
result.te_mf = mTe;

result.SubNet_Full = SubNet_Full;
result.SubNet_List = SubNet_Trimmed;
result.SubNet_MI = SubNet_MI;
result.SubNet_Score = SubNet_PVal;
result.Gene_Name = dataset_info.DatasetTr.Gene_Name;
end

%% -----------------------------  SUB functions ------------------------------------
function valid_subnet = filterSubNetworks_PermExprRndNode(SubNet_MI, pv_tresh, Expression, Lbl, SubNet_Full, GAMMA, n_bin)
n_gene = size(Expression, 2);
Expression = Expression(:, randperm(n_gene));
n_snet = numel(SubNet_MI);
n_rndNet = 100; % According to Author's paper
randNet_MI = zeros(100,1);
fprintf('Filter 1: Generating [%d] random subnetworks using [random seed] on permuted network: ', n_rndNet);
for ri=1:n_rndNet
	showprogress(ri, n_rndNet);
	rnd_i = randi(n_gene, 1, 1);
	[~, randNet_MI(ri)] = trimSubNetwork(SubNet_Full{rnd_i}, Expression, Lbl, GAMMA, n_bin);
end
[phat, ~] = gamfit(randNet_MI, pv_tresh);
[~, pcov] = gamlike(phat, randNet_MI);
valid_subnet = true(n_snet,1);
for si=1:n_snet
    if 1-gamcdf(SubNet_MI(si),phat(1),phat(2),pcov)>pv_tresh
        valid_subnet(si) = 0;
    end
end
fprintf('Filtering finished [%d --> %d].\n', n_snet, sum(valid_subnet));
end

function valid_subnet = filterSubNetworks_PermExprSameNode(SubNet_MI, pv_tresh, Expression, Lbl, SubNet_Full, GAMMA, n_bin)
n_gene = size(Expression, 2);
n_snet = numel(SubNet_MI);
valid_subnet = true(n_snet,1);
n_epoch = 100;
randNet_MI = zeros(n_epoch, n_snet);

fprintf('Filter 2: Generating [%d] random subnetworks per gene using [same gene] as seed on permuted expressions: ', n_epoch);
for ei=1:n_epoch
	showprogress(ei, n_epoch);
    Expression = Expression(:, randperm(n_gene));
    for si=1:n_snet
		[~, randNet_MI(ei, si)] = trimSubNetwork(SubNet_Full{si}, Expression, Lbl, GAMMA, n_bin);
    end
end

fprintf('Filter 2: Testing significance of [%d] sub-networks: ', n_snet);
for si=1:n_snet
	showprogress(si, n_snet);
    [phat, ~] = gamfit(randNet_MI(:, si), pv_tresh);
    [~, pcov] = gamlike(phat, randNet_MI(:, si));
    if 1-gamcdf(SubNet_MI(si),phat(1),phat(2),pcov)>pv_tresh
        valid_subnet(si) = 0;
    end
end
fprintf('Filtering finished [%d --> %d].\n', n_snet, sum(valid_subnet));
end

function [valid_subnet, SubNet_PVal] = filterSubNetworks_PermLabels(SubNet_Trimmed, SubNet_MI, pv_tresh, Expression, Lbl, n_bin)
n_sample = size(Expression, 1);
n_snet = numel(SubNet_MI);
n_epoch = 1000;
valid_subnet = true(n_snet,1);
SubNet_PVal = inf(n_snet, 1);
fprintf('Filter 3: Evaluating subnetworks using permuted labels: ');
for si=1:n_snet
	showprogress(si, n_snet);
	n_node = numel(SubNet_Trimmed{si});
	meta = sum(Expression(:,SubNet_Trimmed{si})/sqrt(n_node),2);
	
	randNet_MI = zeros(n_epoch,1);
	for ei=1:n_epoch
		randNet_MI(ei) = getMI(Lbl(randperm(n_sample)), meta, n_bin);
	end
	
    [phat, ~] = gamfit(randNet_MI, pv_tresh);
    [~, pcov] = gamlike(phat, randNet_MI);
    SubNet_PVal(si) = 1 - gamcdf(SubNet_MI(si),phat(1),phat(2),pcov);
    if SubNet_PVal(si)>pv_tresh
        valid_subnet(si) = 0;
    end
end
fprintf('Filtering finished [%d --> %d].\n', n_snet, sum(valid_subnet));
end

function res_mi = getMI(TrueLbl, PredLbl, n_bin)
bin_cen = linspace(min(PredLbl), max(PredLbl), n_bin);
bin_width = (bin_cen(2) - bin_cen(1))/2;

DisLbl = zeros(size(PredLbl));
for bi=1:n_bin
	DisLbl(PredLbl>=bin_cen(bi)-bin_width & PredLbl<bin_cen(bi)+bin_width, 1) = bi;
end
res_mi = mi(DisLbl, TrueLbl);
end

function [subnet_trimmed, old_MI] = trimSubNetwork(SubNet_Full, xTr, lTr, GAMMA, n_bin)
old_MI = getMI(lTr, xTr(:,SubNet_Full(1)), n_bin);
subnet_trimmed = SubNet_Full(1);
for ni=2:numel(SubNet_Full)
	meta = sum(xTr(:,SubNet_Full(1:ni))/sqrt(ni),2);
	new_MI = getMI(lTr, meta, n_bin);
	if (new_MI - old_MI)/old_MI < GAMMA
		subnet_trimmed = SubNet_Full(1:ni-1);
		break;
	end
	old_MI = new_MI;
end
if isequal(ni,numel(SubNet_Full)) && (new_MI - old_MI)/old_MI >= GAMMA % if no break has happend
	subnet_trimmed = SubNet_Full;
end
end

function Nei_lst = getNetNeighborsBreadthFirst(Neig_cell, Nei_lst, n_neighbor, Depth_ind)
if numel(Nei_lst)>=n_neighbor
	Nei_lst = Nei_lst(1:n_neighbor);
elseif Depth_ind<numel(Nei_lst)
	Depth_ind = Depth_ind + 1;
	Nei_lst = unique([Nei_lst Neig_cell{Nei_lst(Depth_ind)}], 'stable');
	if numel(Nei_lst)>=n_neighbor
		Nei_lst = Nei_lst(1:n_neighbor);
	else
		Nei_lst = getNetNeighborsBreadthFirst(Neig_cell, Nei_lst, n_neighbor, Depth_ind);
	end
	Nei_lst = Nei_lst(1:min([numel(Nei_lst) n_neighbor]));
end
end

