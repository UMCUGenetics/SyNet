function S07_GetPerformance_of_StdPairs()

%% Initialization
cv_info = load('/Users/amin/Technical/My_Works/Deep_Learning/113_Organize_SPADE_Codes/Analysis/08_Perform_NOPs_Restricted_Network/CV_Files/CV_SyNet-SyNet_CVT01_Si01-Ri001.mat');
iTr = cv_info.tr_info.CVInd;
iTe = cv_info.te_info.CVInd;
% [~, sid] = sort(abs(Top_Pair(:,10)), 'Descend');
% Pair_ind = Top_Pair(sid,:);
% Pair_ind = Top_Pair(1,1:2);
Pair_ind = Top_Pair(abs(Top_Pair(:,10))>0.05,:);
Pair_ind = Top_Pair(abs(Top_Pair(:,6))>0.6,:);
gind = unique(Pair_ind(1:10,1:2)', 'Stable');
xTr = ge_data.Gene_Expression(iTr, gind);
xTe = ge_data.Gene_Expression(iTe, gind);
lTr = double(ge_data.Patient_Label(iTr)==1)*2-1;
lTe = double(ge_data.Patient_Label(iTe)==1)*2-1;
dataset_info.DatasetTr.iCvPar = cv_info.tr_info.iCvPar;
MAX_N_SUBNET = 500;
opt_info.lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};
n_gene = size(xTr, 2);


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
fprintf('Training the final model over [%d] features...\n', n_gene);
result = LassoWithCV(@lassoEx, zTr, lTr, zTe, lTe, dataset_info.DatasetTr.iCvPar, opt_info.lasso_opt);
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);

