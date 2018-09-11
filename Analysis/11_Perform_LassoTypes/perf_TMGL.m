function result = perf_TMGL(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
Net_Adj = dataset_info.DatasetTr.Net_Adj;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
Lambda_lst = opt_info.lam_list;
n_lam = size(Lambda_lst, 1);
NeigSize_lst = opt_info.NeigSize_lst;
IndexBestNei = 1;

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);

%% Filter network
fprintf('[i] Limiting genes in the network to [%d] ...\n', opt_info.MAX_N_Gene);
[Net_Adj, ~] = SelectTopFromNet(Net_Adj, 'G', opt_info.MAX_N_Gene);
del_ind = sum(Net_Adj~=0)==0;
Net_Adj(del_ind, :) = [];
Net_Adj(:, del_ind) = [];
zTr(:, del_ind) = [];
zTe(:, del_ind) = [];
Gene_Name = dataset_info.DatasetTr.Gene_Name(~del_ind);

%% Adjast the network
if isfield(opt_info, 'AdjustNet') && opt_info.AdjustNet==1
    fprintf('\n[i] Adjusting network according to TM features.\n');
    n_pair = 7088;
    dsn_fname = sprintf('../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_%s_AvgSynACr.mat', DatasetTr.Net_info.net_source);
    DSN_Info = load(dsn_fname, 'Gene_Name');
    res_path = '../19_Compute_Topological_Measures_For_Network/Saved_Pred/';
    res_ptr = sprintf('PredTM_CL-Lasso_SM-OneGRND_NS-%d_NF-50_ID-*.mat', n_pair);
    res_lst = dir([res_path res_ptr]);
    n_res = numel(res_lst);
    fprintf('[%d] TM prediction results are found.\n', n_res);
    TM_PMat = zeros(n_pair, n_res);
    for ri=1:n_res
        res_name = [res_path res_lst(ri).name];
        %fprintf('Loading TM prediction [%s]\n', res_name);
        res_info = load(res_name);
        if ri==1
            TM_PInfo = res_info.Pair_Info;
        else
            if ~isequal(TM_PInfo, res_info.Pair_Info), error(); end
        end
        TM_PMat(:, ri) = res_info.TM_PLabel;
    end
    TM_PMat = TM_PMat ./ size(TM_PMat,2);
    
    Net_GMap = containers.Map(Gene_Name, 1:numel(Gene_Name));
    MAX_NET_WEIGHT = max(Net_Adj(:));
    TM_N_Added = 0;
    TM_N_Removed = 0;
    for pi=1:n_pair
        if Net_GMap.isKey(DSN_Info.Gene_Name{TM_PInfo(pi,1)}) && Net_GMap.isKey(DSN_Info.Gene_Name{TM_PInfo(pi,2)})
            src_ind = Net_GMap(DSN_Info.Gene_Name{TM_PInfo(pi,1)});
            tar_ind = Net_GMap(DSN_Info.Gene_Name{TM_PInfo(pi,2)});
            if Net_Adj(src_ind, tar_ind)==0 && sum(TM_PMat(pi, :))>=0.75
                Net_Adj(src_ind, tar_ind) = MAX_NET_WEIGHT;
                TM_N_Added = TM_N_Added + 1;
            elseif Net_Adj(src_ind, tar_ind)~=0 && sum(TM_PMat(pi, :))<=-0.75
                Net_Adj(src_ind, tar_ind) = 0;
                TM_N_Removed = TM_N_Removed + 1;
            end
            Net_Adj(tar_ind, src_ind) = Net_Adj(src_ind, tar_ind);
        end
    end
    fprintf('Adjustment finished. Modifications: [%d] added, [%d] removed.\n\n', TM_N_Added, TM_N_Removed);
end

%% Generate Neighbor Sets
fprintf('Generating neighbor sets and subnetworks: \n');
Neig_cell = getNeighborsFromAdj(Net_Adj);

%% Prepare final dataset
NeigIndex_lst = GenerateNeighborSets(Neig_cell, NeigSize_lst(IndexBestNei));
[cmbTr, Group_IndexTr] = GenerateGL_Dataset(zTr, NeigIndex_lst);
[cmbTe,             ~] = GenerateGL_Dataset(zTe, NeigIndex_lst);

%% Train Group lasso
fprintf('Training group lasso over [%d] groups and [%d] features ...\n', numel(NeigIndex_lst), size(cmbTr,2));
lasso_opt = {'lassoType', 'sgt', 'CV', 5, 'iCvPar', Fold_Index, 'relTol', 5e-2, 'lam_list', Lambda_lst, 'n_lC', 15, 'lC_ratio', 1e-2, ...
    'n_lG', 15, 'lG_ratio', 1e-2, 'group_ind', Group_IndexTr, 'lambdaType', 'no-grid', 'paroptions', statset('UseParallel',false), 'verbose', 0};
[opt_B, opt_fit] = lassoEx(cmbTr, lTr, lasso_opt{:});

%% Collect Subnet scores
n_snet = size(Group_IndexTr, 2);
SubNet_Score = zeros(n_snet, 1);
IndexBestLamb = opt_fit.IndexMinMSE;
for si=1:n_snet
    feat_index = Group_IndexTr(1,si):Group_IndexTr(2,si);
    SubNet_Score(si) = mean(abs(opt_B(feat_index, IndexBestLamb)));
end

%% Showing results
tr_auc_lam = zeros(1, n_lam);
te_auc_lam = zeros(1, n_lam);
for i=1:n_lam
    tr_auc_lam(i) = getAUC(lTr, cmbTr*opt_B(:,i), 50);
    te_auc_lam(i) = getAUC(lTe, cmbTe*opt_B(:,i), 50);
end
fprintf('    Index: '); fprintf('%5d  ', 1:n_lam); fprintf('\n');
fprintf('Train AUC: '); fprintf('%0.3f, ', 1-opt_fit.MSE); fprintf('\n');
fprintf(' Test AUC: '); fprintf('%0.3f, ', te_auc_lam); fprintf('\n');
fprintf('Optimal lam is: [%d]\n', IndexBestLamb);

%% Plot results
%{
plot(tr_auc_lam); hold on; plot(te_auc_lam); plot(1-opt_fit.MSE);
plot(opt_fit.IndexMinMSE,0.7, '*');
legend({'Train' 'Test' 'CV err'});
%}

%% Storing results
result.B = opt_B;
result.fit = opt_fit;
result.OptLambda = Lambda_lst(IndexBestLamb, :);
result.OptNeiSize = NeigSize_lst(IndexBestNei);
result.OptNetSize = opt_info.MAX_N_Gene;
result.tr_auc_lam = tr_auc_lam;
result.te_auc_lam = te_auc_lam;
result.tr_auc = tr_auc_lam(IndexBestLamb);
result.te_auc = te_auc_lam(IndexBestLamb);
result.SubNet_List = NeigIndex_lst;
result.SubNet_Score = SubNet_Score;
result.Gene_Name = Gene_Name;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
end

function [SubNet_Trimmed, Neig_cell] = GenerateNeighborSets(Neig_cell, MAX_N_Neighbor)
n_gene = numel(Neig_cell);

%% Generate Subnetworks
SubNet_Full = cell(n_gene, 1);
for gi=1:n_gene
    showprogress(gi, n_gene);
    SubNet_Full{gi} = getNetNeighborsBreadthFirst(Neig_cell, Neig_cell{gi}, MAX_N_Neighbor, 1);
end
fprintf('[%d] Subnetworks are collected in breath first manner.\n', numel(SubNet_Full));

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
fprintf('Subnetworks are purified to [%d] unique ones.\n', n_snet);
end

function [CmbData, Group_Index] = GenerateGL_Dataset(GE_Data, NeigIndex_lst)
n_sample = size(GE_Data, 1);
n_snet = numel(NeigIndex_lst);
MAX_N_Neighbor = max(cellfun('length', NeigIndex_lst));

%% Construct Composite features
CmbData = zeros(n_sample, n_snet*MAX_N_Neighbor);
Group_Index = ones(3, n_snet);
step = 1;
for si=1:n_snet
    Group_Index(1, si) = step;
    for gi=1:numel(NeigIndex_lst{si})
        CmbData(:, step) = GE_Data(:, NeigIndex_lst{si}(gi));
        step = step + 1;
    end
    Group_Index(2, si) = step - 1;
end
CmbData(:, step:end) = [];
end

