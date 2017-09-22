function result = perf_iTaylor(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
lTe = dataset_info.DatasetTe.Patient_Label;
Net_Adj = dataset_info.DatasetTr.Net_Adj;
[n_TrSample, n_gene] = size(xTr);
n_TeSample = size(xTe, 1);
MIN_SUBNET_SIZE = 5;
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Generate Neighbor Sets
fprintf('Generating [%d] neighbor sets and subnetworks.\n', n_gene);
Neig_cell = getNeighborsFromAdj(Net_Adj);

%% Trimming Subnetworks By Size
net_size = cellfun(@(x) numel(x), Neig_cell);
Neig_cell(net_size<=MIN_SUBNET_SIZE) = [];
n_snet = numel(Neig_cell);
fprintf('Subnetworks are trimmed. [%d] are left.\n', n_snet);

%% Sorting Subnetworks By Correlation Score
fprintf('Calculating the correlation scores.\n');
cls_lst = unique(lTr);
n_cls = numel(cls_lst); 
if n_cls~=2, error('This function is desinged for 2 class problems only.\n'); end
cr_scr = zeros(n_snet, 1);
for ni=1:n_snet
    cr_mat = zeros(n_cls, numel(Neig_cell{ni})-1);
    for ci=1:n_cls
        cr_mat(ci,:) = corr(xTr(lTr==cls_lst(ci), Neig_cell{ni}(1)), xTr(lTr==cls_lst(ci), Neig_cell{ni}(2:end)));
    end
    cr_scr(ni,1) = mean(abs(cr_mat(1,:)-cr_mat(2,:)));
end
[~, cr_sind] = sort(cr_scr, 'descend');
Neig_cell = Neig_cell(cr_sind);
cr_scr = cr_scr(cr_sind);

%% Generating Meta-features from Top SubNetworks
fprintf('Generating Meta-features from %d sub-networks.\n', MAX_N_SUBNET);
mTr = zeros(n_TrSample, MAX_N_SUBNET);
mTe = zeros(n_TeSample, MAX_N_SUBNET);
for ni=1:min([MAX_N_SUBNET n_snet])
    grp_list = Neig_cell{ni};
    mTr(:, ni) = mean(repmat(xTr(:, grp_list(1)), 1, numel(grp_list)-1) - xTr(:, grp_list(2:end)), 2);
	mTe(:, ni) = mean(repmat(xTe(:, grp_list(1)), 1, numel(grp_list)-1) - xTe(:, grp_list(2:end)), 2);
end

%% Normalization
zTr = zscore(mTr);
zTe = zscore(mTe);
n_meta = size(zTr, 2);

%% Traning the final model
fprintf('Training the final model over [%d] meta-features...\n', n_meta);
result = LassoWithCV(@lassoEx, zTr, lTr, [], [], dataset_info.DatasetTr.iCvPar, opt_info.lasso_opt);
opt_B = result.B;
opt_fit = result.fit;
clear result
fprintf('Final training is done. [%d] non-zero features identified.\n', sum(abs(opt_B(:, opt_fit.IndexMinMSE))>0));

%% Evaluating the model
vec_B = opt_B(:, opt_fit.IndexMinMSE);
g = zTr*vec_B;
tr_auc = getAUC(lTr, g, 50);

g = zTe*vec_B;
te_auc = getAUC(lTe, g, 50);

fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', te_auc*100);

result.B = opt_B;
result.fit = opt_fit;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.tr_mf = mTr;
result.te_mf = mTe;
result.SubNet_List = Neig_cell;
result.SubNet_Score = cr_scr;
result.Gene_Name = dataset_info.DatasetTr.Gene_Name;
end
