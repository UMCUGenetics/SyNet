function result = perf_NN(dataset_info, opt_info)

%% Initialization
zTr = dataset_info.DatasetTr.Gene_Expression;
zTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
n_gene = size(zTr, 2);
Gene_Name = dataset_info.DatasetTr.Gene_Name;
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
n_Fold = max(Fold_Index);
if isfield(opt_info, 'GridSearch') && opt_info.GridSearch==1
    opt_info.n_Hidden = [1 2 4 8 16 32];
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
    fprintf('Dataset now has Train: [%d x %d], Test: [%d x %d] samples and genes.\n', size(zTr), size(zTe));
end
SubNet_List = num2cell(1:n_gene)';

%% Traning the model
n_Hidden = numel(opt_info.n_Hidden);
if n_Hidden==1
    opt_H = opt_info.n_Hidden;
else
    fprintf('CV over:\n#Hidden: %s\n', num2str(opt_info.n_Hidden,'%d, '));
    grid_AUC = zeros(n_Hidden, 1);
    for hi=1:n_Hidden
        fprintf('[%02d] H=%3d, AUC: ', hi, opt_info.n_Hidden(hi));
        cv_auc = ones(n_Fold, 1)*0.5;
        for fi=1:n_Fold
            iTr = Fold_Index~=fi;
            iTe = Fold_Index==fi;
            try
                NNModel = NN_Train(zTr(iTr,:), lTr(iTr), opt_info.n_Hidden(hi), struct('Verbose', 0));
                cv_auc(fi) = getModelAUC(NNModel, zTr(iTe,:), lTr(iTe));
                fprintf('%0.0f, ', cv_auc(fi)*100);
            catch
                fprintf('--, ');
            end
        end
        grid_AUC(hi) = mean(cv_auc);
        fprintf('Mean is: %0.2f\n', grid_AUC(hi)*100);
    end
    
    Top_AUC = max(grid_AUC(:));
    opt_rowi = find(grid_AUC==Top_AUC, 1);
    opt_H = opt_info.n_Hidden(opt_rowi);
    fprintf('Best parameters are: #Hidden(%d) = %d, AUC=%0.2f%%\n', opt_rowi, opt_H, Top_AUC*100);
end

%% Training NN
fprintf('Traning final NN over #Hidden = %d ...\n', opt_H);
NNModel = NN_Train(zTr, lTr, opt_H, struct('Verbose', 1));

%% Evaluating the model
fprintf('Evaluating the model ...\n');
tr_auc = 1;
te_auc = getModelAUC(NNModel, zTe, lTe);

%% Ranking features
if n_gene>1000
    fprintf('Too many features [%d], feature scoring is ignored ...\n', n_gene);
    result.SubNet_Score = randn(n_gene, 1);
else
    fprintf('Ranking features using NN over #Hidden = %d ...\n', opt_H);
    feat_score = zeros(n_Fold, n_gene);
    for fi=1:n_Fold
        fprintf('Folds %2d: ', fi);
        iTr = Fold_Index~=fi;
        iTe = Fold_Index==fi;
        F_zTr = zTr(iTr,:);
        F_lTr = lTr(iTr);
        F_zTe = zTr(iTe,:);
        F_lTe = lTr(iTe);
        n_Te = numel(F_lTe);
        NNModel = NN_Train(F_zTr, F_lTr, opt_H, struct('Verbose', 0));
        org_auc = getModelAUC(NNModel, F_zTe, F_lTe);
        for gi=1:n_gene
            showprogress(gi, n_gene);
            F_zTe(:,gi) = randn(n_Te,1);
            try
                prt_auc = getModelAUC(NNModel, F_zTe, F_lTe);
                feat_score(fi,gi) = org_auc - prt_auc;
            catch
            end
            F_zTe(:,gi) = zTr(iTe,gi);
        end
    end
    feat_zscr = zscore(feat_score, 0, 2);
    result.SubNet_FeatImp = feat_score;
    result.SubNet_Score = mean(feat_zscr, 1);
    fprintf('Feature scoring finished ...\n');
end

%% Saving results
%result.NNModel = NNModel;
result.opt_H = opt_H;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.SubNet_List = SubNet_List;
result.Gene_Name = Gene_Name;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
end

function NNModel = NN_Train(xTr, lTr, n_Hidden, opt)
NNModel = feedforwardnet(n_Hidden);
NNModel = configure(NNModel, xTr', lTr');
NNModel.trainParam.showWindow = 0;
if isfield(opt, 'Verbose') && opt.Verbose==1
    NNModel.trainParam.showCommandLine = 1;
end
NNModel = train(NNModel, xTr', lTr');
end

function model_auc = getModelAUC(NNModel, xTe, lTe)
pred_lbl = NNModel(xTe');
model_auc = getAUC(lTe, pred_lbl, 50);
end