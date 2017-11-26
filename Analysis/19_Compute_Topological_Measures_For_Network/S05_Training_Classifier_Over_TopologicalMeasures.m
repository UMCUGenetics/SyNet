clc;
clear;

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
n_pair = 10000;
IS_STRICT_CV = 1;
CLS_Name = 'Lasso';

%% Load TM data
load('./Topological_Data/TMData_NS20000_NF90.mat', 'zTM_Data', 'TM_Label', 'TM_Name');

%% Evaluation of classifier
fprintf('Training the final model ...\n');
n_Fold = 50;
Fold_Index = crossvalind('KFold', TM_Label, n_Fold);
lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};
CV_auc = zeros(n_Fold, 1);
Feat_Imp = zeros(n_Fold, n_feature);
for fi=1:n_Fold
    iTr = Fold_Index~=fi;
    iTe = Fold_Index==fi;
    
    if IS_STRICT_CV
        is_in = IsInTrain(Pair_Info(:,1:2), Pair_Info(iTe,1:2));
        iTr(is_in) = 0;
    end
    zTr = zTM_Data(iTr, :);
    zTe = zTM_Data(iTe, :);
    lTr = TM_Label(iTr);
    lTe = TM_Label(iTe);
    
    %Sample_Weights = (lTr(iTr)==-1)/sum(lTr(iTr)==-1) + (lTr(iTr)==1)/sum(lTr(iTr)==1);
    %RFModel = TreeBagger(100, zTr(iTr,:), lTr(iTr), 'Method', 'classification', 'MinLeafSize', opt_info.MinLeafSize(li), ...
    %    'Weights', Sample_Weights, 'Surrogate','off', 'PredictorSelection', 'curvature', 'OOBPredictorImportance', 'off');
    %[~, pred_prob] = predict(RFModel, zTr(iTe,:));
    %[~, pred_lbl] = max(pred_prob, [], 2);
    % cv_auc(fi) = getAUC(lTr(iTe), pred_lbl, 50);
    %iCvPar = crossvalind('KFold', TM_Label(iTr), 2);
    switch CLS_Name
        case 'Lasso'
            [B, fit] = lassoEx(zTr, lTr, lasso_opt{:});
            Feat_Imp(fi, :) = B(:, 8);
            CV_auc(fi) = getAUC(lTe, zTe*Feat_Imp(fi, :)', 50) * 100;
            n_lam = size(B, 2);
            tr_auc_lam = zeros(1, n_lam);
            te_auc_lam = zeros(1, n_lam);
            for i=1:n_lam
                tr_auc_lam(i) = getAUC(lTr, zTr*B(:,i), 50);
                te_auc_lam(i) = getAUC(lTe, zTe*B(:,i), 50);
            end
            fprintf('    Index: '); fprintf('%5d  ', 1:n_lam); fprintf('\n');
            fprintf('Train AUC: '); fprintf('%0.3f, ', tr_auc_lam); fprintf('\n');
            fprintf(' Test AUC: '); fprintf('%0.3f, ', te_auc_lam); fprintf('\n');
    end
    fprintf('[%d/%d] Test performance is [%0.2f%%] AUC.\n', fi, n_Fold, CV_auc(fi));
end
fprintf('Mean is: %0.2f\n', mean(CV_auc));

%% Normalizing feature importance
zFeat_Imp = zscore(Feat_Imp, 0, 2);
mFeat_Imp = mean(zFeat_Imp);
[~, sind] = sort(mFeat_Imp, 'Descend');

%% Plotting
figure('Position', [100 100 1500 700]);
imagesc(zFeat_Imp(:, sind));
set(gca, 'XTick', 1:n_feature, 'XTickLabel', TM_Name(sind), 'XTickLabelRotation', 45);
disp(table(TM_Name(sind(1:5)), mFeat_Imp(sind(1:5))', 'VariableNames', {'FeatureName' 'AvgScore'}));

%% Functions %%%%%%%%%%%%%
function is_in=IsInTrain(Query_Index, Train_Index)
is_in = any(ismember(Query_Index, Train_Index(:)), 2);
fprintf('[i] STRICT CV: Out of [%d] query samples, [%d] should be filtered.\n', size(Query_Index,1), sum(is_in));
end