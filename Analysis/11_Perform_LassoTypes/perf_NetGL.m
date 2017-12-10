function result = perf_NetGL(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
Net_Adj = dataset_info.DatasetTr.Net_Adj;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
n_fold = max(Fold_Index);
Lambda_lst = opt_info.lam_list;
n_lam = size(Lambda_lst, 1);

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

%% Generate Neighbor Sets
fprintf('Generating neighbor sets and subnetworks: \n');
Neig_cell = getNeighborsFromAdj(Net_Adj);

%% Select number of neighbor
NeigSize_lst = [2 3 5 7 10 20];
n_nei = numel(NeigSize_lst);
Grid_auc = zeros(n_nei, n_lam, n_fold);
for ni=1:numel(NeigSize_lst)
    fprintf('[i] Using neighbor size = %d\n', NeigSize_lst(ni));

    %% Generate Group Lasso dataset
    NeigIndex_lst = GenerateNeighborSets(Neig_cell, NeigSize_lst(ni));
    [cmbTr, Group_Index] = GenerateGL_Dataset(zTr, NeigIndex_lst);
    lasso_opt = {'lassoType', 'sgt', 'CV', [], 'relTol', 5e-2, 'lam_list', Lambda_lst, 'n_lC', 15, 'lC_ratio', 1e-2, ...
        'n_lG', 15, 'lG_ratio', 1e-2, 'group_ind', Group_Index, 'lambdaType', 'no-grid', 'paroptions', statset('UseParallel',false), 'verbose', 0};

    for fi=1:n_fold
        %% Train over folds
        fTr = Fold_Index~=fi;
        fTe = Fold_Index==fi;
        
        fprintf('Fold [%02d/%02d]: Training over [#Tr=%4d, #Te=%3d] and [%d] features, ', fi, n_fold, sum(fTr), sum(fTe), size(cmbTr,2));
        [fold_B] = lassoEx(cmbTr(fTr,:), lTr(fTr), lasso_opt{:});
        fold_pred = cmbTr(fTe,:)*fold_B;
        for li=1:n_lam
            Grid_auc(ni,li,fi) = getAUC(lTr(fTe), fold_pred(:,li), 50);
        end
        fprintf('Max AUC is: %0.1f%%\n', max(Grid_auc(ni,:,fi))*100);
    end
    [nei_OAuc, nei_OLam] = max(mean(Grid_auc(ni,:,:),3));
    fprintf('Best AUC is [%0.1f%%] with [%d]th lambda.\n\n', nei_OAuc*100, nei_OLam);
end

%% Display the grid
iCV_AUC = mean(Grid_auc, 3);
fprintf('        %s\n', sprintf('[%0.2d],  ',1:n_lam));
for ni=1:n_nei
    fprintf('%5d: ', NeigSize_lst(ni));
    fprintf('%4.1f%%, ', iCV_AUC(ni,:)*100);
    fprintf('\n');
end
[IndexBestNei, IndexBestLamb] = find(iCV_AUC==max(iCV_AUC(:)),1);
fprintf('Best dataset has [%d] neighbors using [%d]th lambda.\n', NeigSize_lst(IndexBestNei), IndexBestLamb);

%% Prepare final dataset
NeigIndex_lst = GenerateNeighborSets(Neig_cell, NeigSize_lst(IndexBestNei));
[cmbTr, Group_IndexTr] = GenerateGL_Dataset(zTr, NeigIndex_lst);
[cmbTe,             ~] = GenerateGL_Dataset(zTe, NeigIndex_lst);
    
%% Train Group lasso
fprintf('Training group lasso over [%d] groups and [%d] features ...\n', numel(NeigIndex_lst), size(cmbTr,2));
lasso_opt = {'lassoType', 'sgt', 'CV', 5, 'iCvPar', Fold_Index, 'relTol', 5e-2, 'lam_list', Lambda_lst, 'n_lC', 15, 'lC_ratio', 1e-2, ...
	'n_lG', 15, 'lG_ratio', 1e-2, 'group_ind', Group_IndexTr, 'lambdaType', 'no-grid', 'paroptions', statset('UseParallel',false), 'verbose', 1};
[opt_B, opt_fit] = lassoEx(cmbTr, lTr, lasso_opt{:});

%% Collect Subnet scores
SubNet_Score = zeros(n_snet, 1);
for si=1:n_snet
	SubNet_Score(si) = mean(abs(opt_B(Group_IndexTr(1:2,si), IndexBestLamb)));
end

%% Showing results
n_lam = size(opt_B, 2);
tr_auc_lam = zeros(1, n_lam);
te_auc_lam = zeros(1, n_lam);
for i=1:n_lam
	tr_auc_lam(i) = getAUC(lTr, cmbTr*opt_B(:,i), 50);
	te_auc_lam(i) = getAUC(lTe, cmbTe*opt_B(:,i), 50);
end
fprintf('    Index: '); fprintf('%5d  ', 1:n_lam); fprintf('\n');
fprintf('Train AUC: '); fprintf('%0.3f, ', 1-opt_fit.MSE); fprintf('\n');
fprintf(' Test AUC: '); fprintf('%0.3f, ', te_auc_lam); fprintf('\n');
fprintf('Optimal lam is: [%d]\n', opt_fit.IndexMinMSE);
fprintf('Utilized lam is: [%d]\n', IndexBestLamb);

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
