function result = perf_KNN(dataset_info, opt_info)

%% Initialization
zTr = dataset_info.DatasetTr.Gene_Expression;
zTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
n_gene = size(zTr, 2);
Gene_Name = dataset_info.DatasetTr.Gene_Name;
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Sort features using t-test
if isfield(opt_info, 'UseTTest')
    fprintf('Selecting top [%d] features from [%d] genes using t-test.\n', MAX_N_SUBNET, n_gene);
    
    %% Calculating pvalue for genes
    pv_mat = zeros(n_gene, 1);
    for gi=1:n_gene
        [~, pv_mat(gi)] = ttest(zTr(:,gi), lTr);
    end
    
    %% Select top features
    [SubNet_Score, scr_ind] = sort(-log10(pv_mat), 'Descend');
    zTr = zTr(:, scr_ind(1:MAX_N_SUBNET));
    zTe = zTe(:, scr_ind(1:MAX_N_SUBNET));
    Gene_Name = Gene_Name(scr_ind(1:MAX_N_SUBNET));
    n_gene = size(zTr, 2);
    fprintf('Dataset now has Train: [%d x %d], Test: [%d x %d] samples and genes.\n', size(zTr), size(zTe));
else
    SubNet_Score = ones(n_gene, 1);
end
SubNet_List = num2cell(1:n_gene)';

%% Finding the optimal K
if opt_info.K==0
    fprintf('Training the adaptive KNN over [%d] features...\n', n_gene);
    [~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
    [pred, opt_K] = knn(zTr, lTr, zTe, setfield(opt_info, 'iCvPar', Fold_Index));
else
    fprintf('Training KNN using [%d] neighbors over [%d] features ...\n', opt_info.K, n_gene);
    [pred, opt_K] = knn(zTr, lTr, zTe, struct('K', opt_info.K));
end
fprintf('Best K is: %d\n', opt_K);
result.opt_K = opt_K;

%% Evaluating the model
tr_auc = 1;
te_auc = getAUC(lTe, pred, 50);

%% Saving results
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.SubNet_List = SubNet_List;
result.SubNet_Score = SubNet_Score;
result.Gene_Name = Gene_Name;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
end


