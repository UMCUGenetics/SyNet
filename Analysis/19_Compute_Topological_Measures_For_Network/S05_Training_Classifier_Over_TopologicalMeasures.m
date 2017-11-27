clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};
IS_STRICT_CV = 1;
CLS_Name = 'Lasso';
CLS_Name = 'RandomForest';
CLS_Name = 'Regress';
Regex_lst = {
    '^HPRD'
    '^I2D'
    '^STRING'
    '^HBEpith'
    '^(?!HBGland).'
    '^HBGland'
    '.'
    };
n_regex = numel(Regex_lst);
n_Fold = 50;

%% Load TM data
load('./Topological_Data/TMData_NS20000_NF90.mat', 'zTM_Data', 'TM_Label', 'TM_Name', 'Pair_Info');
%qTM_Data = quantilenorm(zTM_Data); % , 'Display', true

%% Evaluation of classifier
Grp_AUC = zeros(n_Fold, n_regex);
Grp_Name = cell(n_regex,1);
for ri=1:n_regex
    In_Grp = cellfun('length', regexp(TM_Name, Regex_lst{ri}, 'once'))==1;
    zData = zTM_Data(:, In_Grp);
    %zData = qTM_Data(:, In_Grp);
    Feat_Name = TM_Name(In_Grp);
    
    %if ri==n_regex
    %    %% Reducing dimension Using PCA
    %    [coeff,score,latent,tsquared,explained,mu] = pca(zData, 'NumComponents', 20);
    %    zData = zData * coeff;
    %    TM_Name = repmat({'PC'},20,1);
    %end
    n_feature = size(zData, 2);
    
    %% Traning the classifier
    fprintf('Training the [%s] model over [%s] ...\n', CLS_Name, Regex_lst{ri});
    Fold_Index = crossvalind('KFold', TM_Label, n_Fold);
    Fold_auc = zeros(n_Fold, 1);
    Feat_Imp = zeros(n_Fold, n_feature);
    for fi=1:n_Fold
        iTr = Fold_Index~=fi;
        iTe = Fold_Index==fi;
        
        if IS_STRICT_CV
            Test_Index = Pair_Info(iTe,1:2);
            is_in = any(ismember(Pair_Info(:,1:2), Test_Index(:)), 2);
            iTr(is_in) = 0;
        end
        zTr = zData(iTr, :);
        zTe = zData(iTe, :);
        lTr = TM_Label(iTr);
        lTe = TM_Label(iTe);
        
        switch CLS_Name
            case 'Regress'
                Feat_Imp(fi, :) = regress(lTr, zTr);
                Fold_auc(fi) = getAUC(lTe, zTe*Feat_Imp(fi, :)', 50) * 100;
            case 'Lasso'
                [B, fit] = lassoEx(zTr, lTr, lasso_opt{:});
                Feat_Imp(fi, :) = B(:, 8);
                Fold_auc(fi) = getAUC(lTe, zTe*Feat_Imp(fi, :)', 50) * 100;
                %{
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
                %}
            case 'RandomForest'
                Sample_Weights = (lTr==-1)/sum(lTr==-1) + (lTr==1)/sum(lTr==1);
                RFModel = TreeBagger(50, zTr, lTr, 'Method', 'classification', 'MinLeafSize', 1, 'NumPrint', 0, ...
                    'Weights', Sample_Weights, 'Surrogate','off', 'PredictorSelection', 'curvature', 'OOBPredictorImportance', 'on');
                [~, pred_prob] = predict(RFModel, zTe);
                [~, pred_lbl] = max(pred_prob, [], 2);
                Fold_auc(fi) = getAUC(lTe, pred_lbl, 50)*100;
                Feat_Imp(fi, :) = RFModel.OOBPermutedPredictorDeltaError;
        end
        fprintf('[%d/%d] Test performance is [%0.2f%%] AUC.\n', fi, n_Fold, Fold_auc(fi));
    end
    fprintf('Mean is: %0.2f\n', mean(Fold_auc));
    
    Grp_AUC(:, ri) = Fold_auc;
    Grp_Name{ri} = sprintf('%s', Regex_lst{ri}(2:end));
end
Grp_Name([5 end]) = {'All but HBGland' 'All networks'};

%% Plotting AUC per group of features
figure('Position', [100 100 700 500]);
hold on
clr_map = hsv(n_regex)*0.8;
for ri=1:n_regex
    box_h = BoxPlotEx(Grp_AUC(:, ri), 'Positions', ri, 'Color', clr_map(ri,:), 'Symbol', '', 'Widths', 0.8);
    set(box_h, 'LineWidth', 1.5);
end
set(gca, 'XTick', 1:n_regex, 'XTickLabel', Grp_Name, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', ...
    'XLim', [0 n_regex+1], 'YLim', [55 95]);
ylabel(sprintf('AUC (across %d folds)', n_Fold));
title('Prediction of SyNet links', 'FontSize', 12);
if IS_STRICT_CV
    output_name = sprintf('./Plots/S05_ClassifierPerformance_%s_Over_TM_NF%d_UseStrictCV.pdf', CLS_Name, n_Fold);
else
    output_name = sprintf('./Plots/S05_ClassifierPerformance_%s_Over_TM_NF%d.pdf', CLS_Name, n_Fold);
end
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [5 4], 'PaperPosition', [0 0 5 4]);
print('-dpdf', '-r300', output_name);

return

%% Normalizing feature importance
zFeat_Imp = zscore(Feat_Imp, 0, 2);
mFeat_Imp = mean(zFeat_Imp);
[~, sind] = sort(mFeat_Imp, 'Descend');

%% Plotting
figure('Position', [100 100 1500 700]);
imagesc(zFeat_Imp(:, sind));
set(gca, 'XTick', 1:n_feature, 'XTickLabel', TM_Name(sind), 'XTickLabelRotation', 45);
disp(table(TM_Name(sind(1:5)), mFeat_Imp(sind(1:5))', 'VariableNames', {'FeatureName' 'AvgScore'}));

