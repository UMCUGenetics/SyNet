function result = perf_SVM(dataset_info, opt_info)

%% Initialization
zTr = dataset_info.DatasetTr.Gene_Expression;
zTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
n_gene = size(zTr, 2);
Gene_Name = dataset_info.DatasetTr.Gene_Name;
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
if isfield(opt_info, 'GridSearch') && opt_info.GridSearch==1
    opt_info.C = logspace(0, 10, 7);
    opt_info.gamma = logspace(0, 10, 7);
end

%% Sort features using t-test
if isfield(opt_info, 'UseTTest')
    fprintf('Selecting top [%d] features from [%d] genes using t-test.\n', MAX_N_SUBNET, n_gene);
    
    %% Calculating pvalue for genes
    pv_mat = zeros(n_gene, 1);
    for gi=1:n_gene
        [~, pv_mat(gi)] = ttest(zTr(:,gi), lTr);
    end
    
    %% Select top features
    [~, scr_ind] = sort(-log10(pv_mat), 'Descend');
    zTr = zTr(:, scr_ind(1:MAX_N_SUBNET));
    zTe = zTe(:, scr_ind(1:MAX_N_SUBNET));
    Gene_Name = Gene_Name(scr_ind(1:MAX_N_SUBNET));
    n_gene = size(zTr, 2);
    fprintf('Dataset now has Train: [%d x %d], Test: [%d x %d] samples and genes.\n', size(zTr), size(zTe));
end
SubNet_List = num2cell(1:n_gene)';

%% Traning the model
n_C = numel(opt_info.C);
n_G = numel(opt_info.gamma);
if n_C==1 && n_G==1
    opt_C = opt_info.C;
    opt_gamma = opt_info.gamma;
else
    n_Fold = max(Fold_Index);
    grid_AUC = zeros(n_C, n_G);
    for ci=1:n_C
        for gi=1:n_G
            fprintf('[%02d/%02d] C=%5.0e, g=%5.0e, AUC: ', ci, gi, opt_info.C(ci), opt_info.gamma(gi));
            cv_auc = ones(n_Fold, 1)*0.5;
            for fi=1:n_Fold
                iTr = Fold_Index~=fi;
                iTe = Fold_Index==fi;
                try
                    SVMModel = SVM(zTr(iTr,:), lTr(iTr), opt_info.kernel, opt_info.C(ci), opt_info.gamma(gi));
                    [~, pred_prob] = predict(SVMModel, zTr(iTe,:));
                    [~, pred_lbl] = max(pred_prob, [], 2);
                    cv_auc(fi) = getAUC(lTr(iTe), pred_lbl, 50);
                    fprintf('%0.0f, ', cv_auc(fi)*100);
                catch
                    fprintf('--, ');
                end
            end
            grid_AUC(ci, gi) = mean(cv_auc);
            fprintf('Mean is: %0.2f\n', grid_AUC(ci,gi)*100);
        end
    end
    
    fprintf(repmat(' ',1,10));
    fprintf('%5.0e   ',  opt_info.gamma);
    fprintf('\n');
    for ci=1:n_C
        fprintf('%5.0e --> ', opt_info.C(ci));
        for gi=1:n_G
            fprintf('%5.1f   ', grid_AUC(ci,gi)*100);
        end
        fprintf('\n');
    end
    
    Top_AUC = max(grid_AUC(:));
    [opt_rowi, opt_coli] = find(grid_AUC==Top_AUC, 1);
    opt_C = opt_info.C(opt_rowi);
    opt_gamma = opt_info.gamma(opt_coli);
    fprintf('Best parameters are: C(%d)=%5.0e, g(%d)=%5.0e, AUC=%0.2f%%\n', opt_rowi, opt_C, opt_coli, opt_gamma, Top_AUC*100);
end

%% Training SVM
fprintf('Traning final [%s] SVM over C=%5.0e and gamma=%5.0e ...\n', opt_info.kernel, opt_C, opt_gamma);
SVMModel = SVM(zTr, lTr, opt_info.kernel, opt_C, opt_gamma);

%% Evaluating the model
fprintf('Evaluating the model ...\n');
tr_auc = 1;
te_auc = getModelAUC(SVMModel, zTe, lTe);

%% Ranking features
fprintf('Ranking features using [%s] SVM over C=%5.0e and gamma=%5.0e ...\nFolds: ', opt_info.kernel, opt_C, opt_gamma);
feat_score = zeros(n_Fold, n_gene);
for fi=1:n_Fold
    fprintf('%d, ', fi);
    iTr = Fold_Index~=fi;
    iTe = Fold_Index==fi;
    F_zTr = zTr(iTr,:);
    F_lTr = lTr(iTr);
    F_lTe = lTr(iTe);
    n_Te = numel(F_lTe);
    SVMModel = SVM(F_zTr, F_lTr, opt_info.kernel, opt_C, opt_gamma);
    org_auc = getModelAUC(SVMModel, zTr(iTe,:), F_lTe);
    for gi=1:n_gene
        try
            F_zTe = zTr(iTe,:);
            F_zTe(:,gi) = randn(n_Te,1);
            prt_auc = getModelAUC(SVMModel, F_zTe, F_lTe);
            feat_score(fi,gi) = org_auc - prt_auc;
        catch
        end
    end
end
feat_zscr = zscore(feat_score')';
fs_average = mean(feat_zscr, 1);
fprintf('\nFeature scoring finished ...\n');

%% Saving results
%result.SVMModel = SVMModel;
result.opt_C = opt_C;
result.opt_gamma = opt_gamma;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.SubNet_List = SubNet_List;
result.SubNet_Score = fs_average;
result.Gene_Name = Gene_Name;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
end

function SVMModel = SVM(X, Y, kernel_name, C, gamma)
if strcmpi(kernel_name, 'rbf')
    SVMModel = fitcsvm(X, Y, 'Verbose', 0, 'KernelFunction', kernel_name, 'BoxConstraint', C, 'KernelScale', gamma);
else
    SVMModel = fitcsvm(X, Y, 'Verbose', 0, 'KernelFunction', kernel_name, 'BoxConstraint', C);
end
end

function model_auc = getModelAUC(SVMModel, xTe, lTe)
[~, pred_prob] = predict(SVMModel, xTe);
[~, pred_lbl] = max(pred_prob, [], 2);
model_auc = getAUC(lTe, pred_lbl, 50);
end