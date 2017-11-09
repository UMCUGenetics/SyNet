function result = perf_RandomForest(dataset_info, opt_info)

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
    opt_info.MinLeafSize = floor(logspace(log10(1), log10(100), 7));
end

%% Sort features using t-test
if isfield(opt_info, 'UseTTest')
    fprintf('Selecting top [%d] features from [%d] genes using t-test.\n', MAX_N_SUBNET, n_gene);
    
    %% Calculating pvalue for genes
    pv_vec = ttest2Ex(zTr, lTr);
    
    %% Select top features
    [~, scr_ind] = sort(-log10(pv_vec), 'Descend');
    zTr = zTr(:, scr_ind(1:MAX_N_SUBNET));
    zTe = zTe(:, scr_ind(1:MAX_N_SUBNET));
    Gene_Name = Gene_Name(scr_ind(1:MAX_N_SUBNET));
    n_gene = size(zTr, 2);
    result.TTestPval = pv_mat(scr_ind);
    fprintf('Dataset now has Train: [%d x %d], Test: [%d x %d] samples and genes.\n', size(zTr), size(zTe));
end
SubNet_List = num2cell(1:n_gene)';

%% Traning the model
n_LeafSize = numel(opt_info.MinLeafSize);
if n_LeafSize==1
    opt_LeafSize = opt_info.MinLeafSize;
else
    n_Fold = max(Fold_Index);
    grid_AUC = zeros(n_LeafSize, 1);
    for li=1:n_LeafSize
        fprintf('[%02d/%02d] Leaf_size=%5d, AUC: ', li, n_LeafSize, opt_info.MinLeafSize(li));
        cv_auc = ones(n_Fold, 1)*0.5;
        for fi=1:n_Fold
            iTr = Fold_Index~=fi;
            iTe = Fold_Index==fi;
            try
                Sample_Weights = (lTr(iTr)==-1)/sum(lTr(iTr)==-1) + (lTr(iTr)==1)/sum(lTr(iTr)==1);
                RFModel = TreeBagger(100, zTr(iTr,:), lTr(iTr), 'Method', 'classification', 'MinLeafSize', opt_info.MinLeafSize(li), ...
                    'Weights', Sample_Weights, 'Surrogate','off', 'PredictorSelection', 'curvature', 'OOBPredictorImportance', 'off');
                [~, pred_prob] = predict(RFModel, zTr(iTe,:));
                [~, pred_lbl] = max(pred_prob, [], 2);
                cv_auc(fi) = getAUC(lTr(iTe), pred_lbl, 50);
                fprintf('%0.0f, ', cv_auc(fi)*100);
            catch
                fprintf('--, ');
            end
        end
        grid_AUC(li) = mean(cv_auc);
        fprintf('Mean is: %0.2f\n', grid_AUC(li)*100);
    end
    
    Top_AUC = max(grid_AUC(:));
    opt_rowi = find(grid_AUC==Top_AUC, 1);
    opt_LeafSize = opt_info.MinLeafSize(opt_rowi);
    fprintf('Best parameters are: MinLeafSize(%d)=%d, AUC=%0.2f%%\n', opt_rowi, opt_LeafSize, Top_AUC*100);
end

%% Training random forest
fprintf('Traning final Random Forest over Train: [%d x %d], Test: [%d x %d] samples and genes ...\n', size(zTr), size(zTe));
Sample_Weights = (lTr==-1)/sum(lTr==-1) + (lTr==1)/sum(lTr==1);
RFModel = TreeBagger(100, zTr, lTr,'Method', 'classification', 'Surrogate','off', 'NumPrint', 10, ...
    'Weights', Sample_Weights, 'PredictorSelection', 'curvature', 'OOBPredictorImportance', 'on', 'MinLeafSize', opt_LeafSize);

%% Evaluating the model
fprintf('Evaluating the model ...\n');
tr_auc = 1;
[~, pred_prob] = predict(RFModel, zTe);
[~, pred_lbl] = max(pred_prob, [], 2);
te_auc = getAUC(lTe, pred_lbl, 50);

%% Saving results
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.SubNet_List = SubNet_List;
result.SubNet_Score = RFModel.OOBPermutedPredictorDeltaError;
result.Gene_Name = Gene_Name;
result.opt_LeafSize = opt_LeafSize;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
end

