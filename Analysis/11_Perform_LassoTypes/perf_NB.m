function result = perf_NB(dataset_info, opt_info)

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
    pv_vec = ttest2Ex(zTr, lTr);
    
    %% Select top features
    [SubNet_Score, scr_ind] = sort(-log10(pv_vec), 'Descend');
    zTr = zTr(:, scr_ind(1:MAX_N_SUBNET));
    zTe = zTe(:, scr_ind(1:MAX_N_SUBNET));
    Gene_Name = Gene_Name(scr_ind(1:MAX_N_SUBNET));
    n_gene = size(zTr, 2);
    fprintf('Dataset now has Train: [%d x %d], Test: [%d x %d] samples and genes.\n', size(zTr), size(zTe));
else
    SubNet_Score = ones(n_gene, 1);
end
SubNet_List = num2cell(1:n_gene)';

%% Training classifier
fprintf('Training the Naive Bayes over [%d] features...\n', n_gene);
[pred, ~] = nb(zTr, lTr, zTe, opt_info);

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

function [label, Posterior] = nb(xTr, lTr, zTe, nb_opt)
Sample_Weights = (lTr==-1)/sum(lTr==-1) + (lTr==1)/sum(lTr==1);
NB_Mdl = fitcnb(xTr, lTr, 'Distribution', 'normal', 'Weights', Sample_Weights);
[label, Posterior, Cost] = predict(NB_Mdl, zTe);
end
