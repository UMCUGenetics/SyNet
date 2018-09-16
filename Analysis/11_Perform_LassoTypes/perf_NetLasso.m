function result = perf_NetLasso(dataset_info, opt_info)

%% Initialization
n_lam = 20;
lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', n_lam, 'lC_ratio', 1e-2, 'verbose', 0};
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
Net_Adj = dataset_info.DatasetTr.Net_Adj;
[~, ~, Fold_Index] = unique(dataset_info.DatasetTr.iCvPar, 'Stable');
n_iFold = max(Fold_Index);
N_Gene_lst = [50 100 300 500 700 1500 3000];
n_net = numel(N_Gene_lst);

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);
clear xTr xTe

%% Filter network
degree_lst = sum(Net_Adj>0, 1);
if isfield(opt_info, 'Net_MinNEdge')
    fprintf('Minimum edge is requested, filtering network to have at least [%d] edges.\n', opt_info.Net_MinNEdge);
    del_ind = degree_lst < opt_info.Net_MinNEdge;
else
    fprintf('Filtering network to have connected nodes only.\n');
    del_ind = degree_lst == 0;
end
Net_Adj(del_ind, :) = [];
Net_Adj(:, del_ind) = [];
zTr(:, del_ind) = [];
zTe(:, del_ind) = [];
Gene_Name = dataset_info.DatasetTr.Gene_Name(~del_ind);
fprintf('#Genes in the network is [%d].\n', size(Net_Adj, 1));
clear del_ind degree_lst

%% Cross-validation
fprintf('[i] Grid search among [%s] options.\n', num2str(N_Gene_lst, '%d '));
grid_auc = zeros(n_net, n_lam, n_iFold);
for ni=1:n_net
    fprintf('Training over [%d] genes ...\n', N_Gene_lst(ni));
    
    %% Filter for net
    [iNet_Adj, ~] = SelectTopFromNet(Net_Adj, 'G', N_Gene_lst(ni));
    del_ind = sum(iNet_Adj~=0)==0;
    iData = zTr(:, ~del_ind);
    
    %% Train over folds
    for fi=1:n_iFold
        fTr = Fold_Index~=fi;
        fTe = Fold_Index==fi;
        
        fprintf('Fold [%02d/%02d]: Training over [#Tr=%4d, #Te=%3d] and [%d] genes, ', fi, n_iFold, sum(fTr), sum(fTe), size(iData,2));
        [fold_B, ~] = lassoEx(iData(fTr,:), lTr(fTr), lasso_opt{:});
        fold_pred = iData(fTe,:)*fold_B;
        for li=1:n_lam
            grid_auc(ni,li,fi) = getAUC(lTr(fTe), fold_pred(:,li), 50);
        end
        fprintf('Max AUC is: %0.1f%%\n', max(grid_auc(ni,:,fi))*100);
    end
    [net_auc, net_auc_ind] = max(mean(grid_auc(ni,:,:),3));
    fprintf('Best AUC is [%0.1f%%] with [%d]th lambda\n', net_auc*100, net_auc_ind);
end

%% Display the grid
Net_AUC = mean(grid_auc, 3);
fprintf('        %s\n', sprintf('[%0.2d],  ',1:n_lam));
for ni=1:n_net
    fprintf('%5d: ', N_Gene_lst(ni));
    fprintf('%4.1f%%, ', Net_AUC(ni,:)*100);
    fprintf('\n');
end
fprintf('Mean of AUC per Lamb: '); fprintf('%0.1f%%, ', mean(Net_AUC,1)*100); fprintf('\n');
fprintf('Mean of AUC per Nets: '); fprintf('%0.1f%%, ', mean(Net_AUC,2)*100); fprintf('\n');
[IndexBestNet, IndexBestLamb] = find(Net_AUC==max(Net_AUC(:)),1);
% [~, IndexBestLamb] = max(mean(Net_AUC,1));
% [~, IndexBestNet] =  max(mean(Net_AUC,2));
fprintf('[i] Best net has [%d] genes and best lambda is [%d].\n\n', N_Gene_lst(IndexBestNet), IndexBestLamb);

%% Select data
[iNet_Adj, Net_Threshold] = SelectTopFromNet(Net_Adj, 'G', N_Gene_lst(IndexBestNet));
del_ind = sum(iNet_Adj~=0)==0;
Tr_Data = zTr(:, ~del_ind);
Te_Data = zTe(:, ~del_ind);
Gene_Name = Gene_Name(~del_ind);
n_gene = size(Tr_Data, 2);

%% Traning the final model
fprintf('Training the final model over [%d] features...\n', n_gene);
result = LassoWithCV(@lassoEx, Tr_Data, lTr, Te_Data, lTe, dataset_info.DatasetTr.iCvPar, opt_info.lasso_opt);
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);

%% Saving the output
result.UsedTrSamples = dataset_info.DatasetTr.UsedSample;
result.UsedTeSamples = dataset_info.DatasetTe.UsedSample;
result.Net_Threshold = Net_Threshold;
result.BestLambIndex = IndexBestLamb;
result.BestNetwork = N_Gene_lst(IndexBestNet);
result.SubNet_List = num2cell(randperm(n_gene))';
result.SubNet_Score = rand(n_gene,1);
result.Gene_Name = Gene_Name;
end


