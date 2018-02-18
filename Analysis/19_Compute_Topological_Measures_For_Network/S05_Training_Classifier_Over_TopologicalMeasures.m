clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
IS_STRICT_CV = 1;
% CLS_Name = 'Regress';
% CLS_Name = 'Lasso';
CLS_Name = 'RandomForest';
% CLS_Name = 'SVM-RBF';
n_Fold = 50;
% Shuff_Method = 'LnkShuff';
Shuff_Method = 'OneGRND';
n_pair = 7088;
% n_pair = 100000;
% n_pair = 200000;
Regex_lst = {
    '^HumanInt'
    '^BioPlex'
    '^BioGRID'
    '^IntAct'
    '^STRING'
    '^HBOvary'
    '^HBBrain'
    '^HBKidney'
    '^HBLympNode'
%     '^(?!HBGland).'
    '^HBGland'
    'DirectConnection|ShortestPath|Jaccard'
    '^.*'
    };
n_regex = numel(Regex_lst);

%% Load TM data
tmdata_name = sprintf('./Topological_Data/TMData-%s_NS%d_NF210.mat', Shuff_Method, n_pair);
load(tmdata_name, 'TM_Data_z', 'TM_Label', 'TM_Name', 'Pair_Info');
%qTM_Data = quantilenorm(zTM_Data); % , 'Display', true

%% Evaluation of classifier
Grp_AUC = zeros(n_Fold, n_regex);
Grp_Name = cell(n_regex,1);
for ri=1:n_regex
    In_Grp = cellfun('length', regexp(TM_Name, Regex_lst{ri}, 'once'))==1;
    zData = TM_Data_z(:, In_Grp);
    %zData = qTM_Data(:, In_Grp);
    Feat_Name = TM_Name(In_Grp);
    TM_PLabel = zeros(n_pair, 1);
    %if ri==n_regex
    %    %% Reducing dimension Using PCA
    %    [coeff,score,latent,tsquared,explained,mu] = pca(zData, 'NumComponents', 20);
    %    zData = zData * coeff;
    %    TM_Name = repmat({'PC'},20,1);
    %end
    n_feature = size(zData, 2);
    if n_feature<1, error(); end
    
    %% Traning the classifier
    fprintf('Training [%s] over [%15s]: ', CLS_Name, Regex_lst{ri});
    Fold_Index = crossvalind('KFold', TM_Label, n_Fold);
    Fold_auc = zeros(n_Fold, 1);
    n_Sample = zeros(n_Fold, 2);
    Feat_Coeff = zeros(n_Fold, n_feature);
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
        n_Sample(fi, :) = [sum(iTr) sum(iTe)];
        
        switch CLS_Name
            case 'Regress'
                warning off
                B = regress(lTr, zTr);
                warning on
                Fold_auc(fi) = getAUC(lTe, zTe*B, 50) * 100;
                Feat_Coeff(fi, :) = B;
            case 'Lasso'
                lasso_opt = {'lassoType', 't', 'CV', 5, 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};
                [B, fit] = lassoEx(zTr, lTr, lasso_opt{:});
                opt_L = fit.IndexMinMSE;
                Feat_Coeff(fi, :) = B(:, opt_L);
                pred_vec = zTe*B(:,opt_L);
                [te_auc, optThresh] = getAUC(lTe, pred_vec, 50);
                Fold_auc(fi) = te_auc * 100;
                TM_PLabel(iTe) = (pred_vec >= optThresh)*2-1;
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
                Feat_Coeff(fi, :) = RFModel.OOBPermutedPredictorDeltaError;
            case 'SVM-RBF'
                AUC_mat = [];
                C_lst = logspace(-1, 3, 20)';
                G_lst = logspace(1, 3, 20)';
                for ci=1:numel(C_lst)
                    for gi=1:numel(G_lst)
                        SVMModel = fitcsvm(zTr, lTr, 'Verbose', 0, 'KernelFunction', 'rbf', 'BoxConstraint', C_lst(ci), 'KernelScale', G_lst(gi));
                        [~, pred_prob] = predict(SVMModel, zTe);
                        [~, pred_lbl] = max(pred_prob, [], 2);
                        AUC_mat(ci,gi) = getAUC(lTe, pred_lbl, 50)*100;
                    end
                end
            otherwise
                error();
        end
%         fprintf('[%d/%d] Test performance is [%0.2f%%] AUC.\n', fi, n_Fold, Fold_auc(fi));
    end
    fprintf('[%d] folds, Median #Tr=%3.0f, #Te=%3.0f, Median AUC is: %0.2f\n', n_Fold, median(n_Sample), median(Fold_auc));
    
    Grp_AUC(:, ri) = Fold_auc;
    Grp_Name{ri} = sprintf('%s', Regex_lst{ri});
    
    if strcmp(Regex_lst{ri}, '^.*')
        sav_name = sprintf('./Saved_Pred/PredTM_CL-%s_SM-%s_NS-%d_NF-%d_ID-%06.0f.mat', CLS_Name, Shuff_Method, n_pair, n_Fold, rand*1e6);
        fprintf('Saving result of full TM pred in [%s]\n', sav_name);
        save(sav_name, 'TM_PLabel', 'TM_Label', 'Pair_Info');
    end
end
% Grp_Name([5 end]) = {'All but HBGland' 'All networks'};
Grp_Name(end) = {'All networks'};

%% Generate output name
output_name = sprintf('./Plots/S05_PredictTM_CL-%s_SM-%s_NS-%d_NF-%d', CLS_Name, Shuff_Method, n_pair, n_Fold);
if IS_STRICT_CV
    output_name = [output_name '_StrictCV']; 
end

%% Plotting AUC per group of features
close all
figure('Position', [100 100 700 500]);
hold on
clr_map = hsv(n_regex)*0.8;
for ri=1:n_regex
    met_clr = getColor(Grp_Name{ri}(2:end));
    box_h = BoxPlotEx(Grp_AUC(:, ri), 'Positions', ri, 'Color', met_clr, 'Symbol', '', 'Widths', 0.7);
    set(box_h, 'LineWidth', 1.5);
end
ylim([50 100]);
y_tick = get(gca, 'YTick');
y_tick_label = arrayfun(@(y) sprintf('%0.0f%%', y), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_regex, 'XTickLabel', Grp_Name, 'XTickLabelRotation', 20, 'XLim', [0 n_regex+1], 'FontWeight', 'Bold', ...
    'YTick', y_tick, 'YTickLabel', y_tick_label);
ylabel(sprintf('AUC (across %d folds)', n_Fold));
title('Prediction of SyNet links', 'FontSize', 12);

% return
%% Saving AUC plot
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [6 4], 'PaperPosition', [0 0 6 4]);
print('-dpdf', '-r300', [output_name '_AUC.pdf']);

%% Normalizing feature importance
if ismember(CLS_Name, {'RandomForest'})
    Feat_Nrm = Feat_Coeff - min(Feat_Coeff(:));
else
    Feat_Nrm = Feat_Coeff;
end
zFeat_Imp = zscore(abs(Feat_Nrm), 0, 2);
% zFeat_Imp = quantilenorm(Feat_Imp', 'Median', 1)';
% imagesc(zFeat_Imp);
avg_Imp = mean(zFeat_Imp, 1);
[~, sind] = sort(avg_Imp, 'Descend');

%% Plotting Feature importance matrix
close all
figure('Position', [100 100 1500 700]);
% imagesc(Feat_Coeff(:, sind));
imagesc(zFeat_Imp(:, sind(1:100)));
set(gca, 'XTick', 1:n_feature, 'XTickLabel', TM_Name(sind), 'XTickLabelRotation', 45, 'TickLength', [0 0]);
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [30 4], 'PaperPosition', [0 0 30 4]);
print('-dpdf', '-r300', [output_name '_FImpMat.pdf']);

%% Save feature importance table
Feat_Imp_tbl = table(TM_Name(sind), mean(Feat_Coeff(:, sind), 1)', 'VariableNames', {'FeatureName' 'AvgScore'});
out_str = evalc('disp(Feat_Imp_tbl)');
fid = fopen([output_name '_FImpTable.txt'], 'w');
fprintf(fid, out_str);
fclose(fid);

%% Save results
save([output_name '_Matlab.mat'], 'CLS_Name', 'Regex_lst', 'Feat_Coeff', 'Grp_AUC', 'Grp_Name', 'n_Sample', 'n_Fold');
