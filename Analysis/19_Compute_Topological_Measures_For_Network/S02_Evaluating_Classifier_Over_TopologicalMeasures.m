clc;
clear;

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
n_pair = 10000;
IS_STRICT_CV = 1;

%% Load labels
SyNet_info = load('../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat', 'PP_Info', 'NP_Info', 'Gene_Name');
Pair_Info = [
    SyNet_info.PP_Info(1:n_pair,:)  ones(n_pair, 1);
    SyNet_info.NP_Info(1:n_pair,:) -ones(n_pair, 1);
    ];
TM_Label = double(Pair_Info(:,16));

%% Collect data
net_lst = {'STRING', 'HPRD', 'HBEpith', 'HBGland'};
n_net = numel(net_lst);
tm_lst = {'ShortestPath' 'PageRank-FB0.95' 'PageRank-FB0.85' 'PageRank-FB0.75' 'PageRank-FB0.65' 'Eigenvector' 'Degree' 'Closeness' 'Betweenness'};
n_tm = numel(tm_lst);
TM_Data = zeros(n_pair*2, n_net*n_tm*2);
TM_Name = cell(n_net*n_tm*2, 1);
step = 1;
for ni=1:n_net
    for ti=1:n_tm
        data_name = sprintf('./Topological_Results/TM_%s_NP%06d_%s.mat', net_lst{ni}, n_pair, tm_lst{ti});
        fprintf('Loading feature from [%s].\n', data_name);
        Data_Info = load(data_name);
        if ~isequal(SyNet_info.Gene_Name, Data_Info.Ref_GeneName), error(); end
        
        TM_Name{step}   = sprintf('%s-%s-A', net_lst{ni}, tm_lst{ti});
        TM_Data(:,step)   = Data_Info.Pair_AvgScore;
        TM_Name{step+1} = sprintf('%s-%s-D', net_lst{ni}, tm_lst{ti});
        TM_Data(:,step+1) = Data_Info.Pair_DifScore;
        step = step + 2;
    end
end
[n_sample, n_feature] = size(TM_Data);

%% Normalization
%imagesc(TM_Data);
for fi=1:n_feature
    is_inf = TM_Data(:,fi)==inf;
    if sum(is_inf)>n_sample*0.5
        fprintf('[w] Warning: Large number of samples in [%s] are inf [%d/%d] (%0.2f%%).\n', TM_Name{fi}, sum(is_inf), n_sample, sum(is_inf)*100/n_sample);
    end
    TM_Data(is_inf,fi) = max(TM_Data(~is_inf,fi))*1.01;
end
if any(isnan(TM_Data(:))) || any(TM_Data(:)==-inf), error(); end
zTM_Data = zscore(TM_Data);
%imagesc(zTM_Data);

%% Evaluation of classifier
fprintf('Training the final model ...\n');
n_Fold = 5;
Fold_Index = crossvalind('KFold', TM_Label, n_Fold);
lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};
CV_auc = zeros(n_Fold, 1);
Feat_Imp = zeros(n_Fold, n_feature);
for fi=1:n_Fold
    iTr = Fold_Index~=fi;
    iTe = Fold_Index==fi;
    
    %Sample_Weights = (lTr(iTr)==-1)/sum(lTr(iTr)==-1) + (lTr(iTr)==1)/sum(lTr(iTr)==1);
    %RFModel = TreeBagger(100, zTr(iTr,:), lTr(iTr), 'Method', 'classification', 'MinLeafSize', opt_info.MinLeafSize(li), ...
    %    'Weights', Sample_Weights, 'Surrogate','off', 'PredictorSelection', 'curvature', 'OOBPredictorImportance', 'off');
    %[~, pred_prob] = predict(RFModel, zTr(iTe,:));
    %[~, pred_lbl] = max(pred_prob, [], 2);
    % cv_auc(fi) = getAUC(lTr(iTe), pred_lbl, 50);
    iCvPar = crossvalind('KFold', TM_Label(iTr), 5);
    result = LassoWithCV(@lassoEx, zTM_Data(iTr, :), TM_Label(iTr), zTM_Data(iTe, :), TM_Label(iTe), iCvPar, lasso_opt);
    Feat_Imp(fi, :) = result.B(:,result.fit.IndexMinMSE);
    CV_auc(fi) = result.te_auc*100;
    fprintf('[%d/%d] Test performance is [%0.2f%%] AUC.\n', fi, n_Fold, CV_auc(fi));
end
fprintf('Mean is: %0.2f\n', mean(CV_auc));

%% Normalizing feature importance
zFeat_Imp = zscore(Feat_Imp, 0, 2);
imagesc(zFeat_Imp);
mFeat_Imp = mean(zFeat_Imp);
[~, sind] = sort(mFeat_Imp, 'Descend');

table(TM_Name(sind(1:10)), mFeat_Imp(sind(1:10))', 'VariableNames', {'FeatureName' 'AvgScore'})

