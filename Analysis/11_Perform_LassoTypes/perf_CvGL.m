function result = perf_CvGL(dataset_info, opt_info)

%% Initialization
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if isfield(opt_info, 'UseParallel') && opt_info.UseParallel
    delete(gcp('nocreate'));
    parpool('local', opt_info.UseParallel);
    fprintf('Parallel pool is loaded with [%d] workers.\n', poolobj.NumWorkers);
else
    fprintf('Parallel pool is not selected, turning off.');
    ps = parallel.Settings;
    ps.Pool.AutoCreate = false;
end
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
oNet_Adj = dataset_info.DatasetTr.Net_Adj;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
n_fold = max(Fold_Index);
Lambda_lst = opt_info.lam_list;
n_lam = size(Lambda_lst, 1);
NeigSize_lst = [2 3 5 7 10];
n_nei = numel(NeigSize_lst);
NGene_lst = [100 300 700 1000 1500 3000];
n_gene = numel(NGene_lst);

%% Normalization
fprintf('Normalizing data ...\n');
oTr = zscore(xTr);
oTe = zscore(xTe);

%% Inner cross-validation
Total_auc = nan(n_gene, n_nei, n_lam, n_fold);
fprintf('[i] Grid search among:\n\t#Genes [%s]\n\tNeighbor sizes [%s].\n', num2str(NGene_lst, '%d '), num2str(NeigSize_lst, '%d '));
for gi=1:n_gene
    
    %% Filter network
    fprintf('%s\n', repmat('.',1,50));
    [Net_Adj, ~] = SelectTopFromNet(oNet_Adj, 'G', NGene_lst(gi));
    del_ind = sum(Net_Adj~=0)==0;
    zTr = oTr(:, ~del_ind);
    zTe = oTe(:, ~del_ind);
    Net_Adj(del_ind, :) = [];
    Net_Adj(:, del_ind) = [];
    Neig_cell = getNeighborsFromAdj(Net_Adj);
    
    %% Select number of neighbor
    for ni=1:n_nei
        fprintf('\tUsing neighbor size = %d, ', NeigSize_lst(ni));
        
        %% Generate Group Lasso dataset
        NeigIndex_lst = GenerateNeighborSets(Neig_cell, NeigSize_lst(ni), 0);
        [cmbTr, Group_Index] = GenerateGL_Dataset(zTr, NeigIndex_lst);
        lasso_opt = {'lassoType', 'sgt', 'CV', [], 'relTol', 5e-2, 'lam_list', Lambda_lst, 'n_lC', 15, 'lC_ratio', 1e-2, ...
            'n_lG', 15, 'lG_ratio', 1e-2, 'group_ind', Group_Index, 'lambdaType', 'no-grid', 'paroptions', statset('UseParallel',false), 'verbose', 0};
        
        %% Train over folds
        parfor fi=1:n_fold
            %fprintf('I am fold [%d]\n', fi);
            fTr = Fold_Index~=fi;
            fTe = Fold_Index==fi;
            
            fold_B = lassoEx(cmbTr(fTr,:), lTr(fTr), lasso_opt{:});
            fold_pred = cmbTr(fTe,:)*fold_B;
            for li=1:n_lam
                Total_auc(gi,ni,li,fi) = getAUC(lTr(fTe), fold_pred(:,li), 50);
            end
        end
        afl_auc = squeeze(mean(Total_auc(gi,ni,:,:),4));
        [nei_OAuc, nei_OLam] = max(afl_auc);
        fprintf('Best Train AUC is [%0.1f%%] with [%d]th lambda.\n', nei_OAuc*100, nei_OLam);
        
        cmbTe = GenerateGL_Dataset(zTe, NeigIndex_lst);
        fold_B = lassoEx(cmbTr, lTr, lasso_opt{:});
        te_pred = cmbTe*fold_B;
        te_auc = zeros(1, n_lam);
        for li=1:n_lam
            te_auc(1,li) = getAUC(lTe, te_pred(:,li), 50);
        end
        fprintf('   Index: '); fprintf('%5d  ', 1:n_lam); fprintf('\n');
        fprintf('  CV AUC: '); fprintf('%0.1f%%, ', afl_auc*100); fprintf('\n');
        fprintf('Test AUC: '); fprintf('%0.1f%%, ', te_auc*100); fprintf('\n');
    end
end
if any(isnan(Total_auc(:))), warning('Found nan values in AUC mat.'); disp(Total_auc); end

%% Display the grid
Grid_AUC = mean(Total_auc, 4);
for ni=1:n_nei
    fprintf('Neighbor size: [%d]\n', NeigSize_lst(ni));
    fprintf('        %s\n', sprintf('[%0.2d],  ',1:n_lam));
    for gi=1:n_gene
        fprintf('%5d: ', NGene_lst(gi));
        fprintf('%4.1f%%, ', Grid_AUC(gi,ni,:)*100);
        fprintf('\n');
    end
    fprintf('\n');
end
[Opt_AUC, Opt_Index] = max(Grid_AUC(:));
[IndexBestNGene, IndexBestNei, IndexBestLamb] = ind2sub(size(Grid_AUC), Opt_Index);
fprintf('Best dataset has [%0.1f%%] AUC with:\n\t[%d] genes,\n\t[%d] neighbors,\n\t[%d]th lambda.\n', Opt_AUC*100, NGene_lst(IndexBestNGene), NeigSize_lst(IndexBestNei), IndexBestLamb);

%% Prepare final dataset
fprintf('[i] Limiting genes in the network to [%d] ...\n', NGene_lst(IndexBestNGene));
[Net_Adj, ~] = SelectTopFromNet(oNet_Adj, 'G', NGene_lst(IndexBestNGene));
del_ind = sum(Net_Adj~=0)==0;
Net_Adj(del_ind, :) = [];
Net_Adj(:, del_ind) = [];
zTr = oTr(:, ~del_ind);
zTe = oTe(:, ~del_ind);
Gene_Name = dataset_info.DatasetTr.Gene_Name(~del_ind);

%% Generate Neighbor Sets
fprintf('Generating neighbor sets with [%d] size ... \n', NeigSize_lst(IndexBestNei));
Neig_cell = getNeighborsFromAdj(Net_Adj);
NeigIndex_lst = GenerateNeighborSets(Neig_cell, NeigSize_lst(IndexBestNei));
[cmbTr, Group_IndexTr] = GenerateGL_Dataset(zTr, NeigIndex_lst);
[cmbTe,             ~] = GenerateGL_Dataset(zTe, NeigIndex_lst);

%% Train Group lasso
fprintf('Training group lasso over [%d] groups and [%d] features ...\n', numel(NeigIndex_lst), size(cmbTr,2));
lasso_opt = {'lassoType', 'sgt', 'CV', 5, 'iCvPar', Fold_Index, 'relTol', 5e-2, 'lam_list', Lambda_lst, 'n_lC', 15, 'lC_ratio', 1e-2, ...
    'n_lG', 15, 'lG_ratio', 1e-2, 'group_ind', Group_IndexTr, 'lambdaType', 'no-grid', 'paroptions', statset('UseParallel',false), 'verbose', 1};
[opt_B, opt_fit] = lassoEx(cmbTr, lTr, lasso_opt{:});

%% Collect Subnet scores
n_snet = size(Group_IndexTr, 2);
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
plot(opt_fit.IndexMinMSE, 1-opt_fit.MSE(opt_fit.IndexMinMSE), '*');
legend({'Train' 'Test' 'CV err'});
%}

%% Storing results
result.B = opt_B;
result.fit = opt_fit;
result.OptLambda = Lambda_lst(IndexBestLamb, :);
result.OptNeiSize = NeigSize_lst(IndexBestNei);
result.OptNetSize = NGene_lst(IndexBestNGene);
result.tr_auc_lam = tr_auc_lam;
result.te_auc_lam = te_auc_lam;
result.tr_auc = tr_auc_lam(IndexBestLamb);
result.te_auc = te_auc_lam(IndexBestLamb);
result.SubNet_List = NeigIndex_lst;
result.SubNet_Score = SubNet_Score;
result.Gene_Name = Gene_Name;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
if isfield(opt_info, 'UseParallel') && opt_info.UseParallel
    fprintf('Shutting down the parallel pool.\n');
    delete(gcp('nocreate'));
end
end

function [SubNet_Trimmed, Neig_cell] = GenerateNeighborSets(Neig_cell, MAX_N_Neighbor, Verbose)
if ~exist('Verbose', 'var'), Verbose=1; end
n_gene = numel(Neig_cell);

%% Generate Subnetworks
SubNet_Full = cell(n_gene, 1);
for gi=1:n_gene
    if Verbose
        showprogress(gi, n_gene);
    end
    SubNet_Full{gi} = getNetNeighborsBreadthFirst(Neig_cell, Neig_cell{gi}, MAX_N_Neighbor, 1);
end
if Verbose
    fprintf('[%d] Subnetworks are collected in breath first manner.\n', numel(SubNet_Full));
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
if Verbose
    fprintf('Subnetworks are purified to [%d] unique ones.\n', n_snet);
end
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
