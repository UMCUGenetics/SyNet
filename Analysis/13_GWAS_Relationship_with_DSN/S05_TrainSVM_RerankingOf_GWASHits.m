clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/getAUC/');
% addpath('../../../../Useful_Sample_Codes/OScore/');
addpath('../_Utilities/');
n_fold = 5;
DEBUG = 0;

%% Loading gene names from GE data
GEPath = getPath('SyNet');
fprintf('Loading gene names from expression data: [%s]\n', GEPath);
GE_Info = load(GEPath, 'Gene_Expression', 'Gene_Name');
if DEBUG
    GE_Info.Gene_Expression = GE_Info.Gene_Expression(:, 1:100);
    GE_Info.Gene_Name = GE_Info.Gene_Name(1:100);
end
n_gene = numel(GE_Info.Gene_Name);

%% Load GWAS hits over Cohort
fid = fopen('./DSN_iCOGS_Hits/iCOGS_Hits_Genes_MD10.0k.tsv', 'r');
% Id	-Log10(pval)	#Hit	#Hit/Size
GWAS_Info = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
fclose(fid);
if DEBUG
    GWAS_Info{1} = GWAS_Info{1}(1:100);
    GWAS_Info{2} = GWAS_Info{2}(1:100);
    GWAS_Info{3} = GWAS_Info{3}(1:100);
    GWAS_Info{4} = GWAS_Info{4}(1:100);
end
if ~isequal(GE_Info.Gene_Name, GWAS_Info{1}), error(); end

%% Loading cancer genes
fid = fopen('./Census_data/Census.tsv', 'r');
% Gene Symbol, Tier, Tumour Types(Somatic), Role in Cancer 
Census_info = textscan(fid, '%s%s%s%s', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
if ~feof(fid), error(); end
fclose(fid);
is_in = ismember(Census_info{1}, GE_Info.Gene_Name);
Census_info = [Census_info{1}(is_in) Census_info{2}(is_in) Census_info{3}(is_in) Census_info{4}(is_in)];
Census_isCancer = double(ismember(GE_Info.Gene_Name, Census_info(:,1)));

%% Loading SyNet
dsn_name = ['../01_Pairwise_Evaluation_of_Genes/Network_Files/DSN_SyNet.mat'];
fprintf('Loading [%s] network.\n', dsn_name);
DSN_info = load(dsn_name, 'Pair_AUC', 'Gene_Name');
if DEBUG
    DSN_info.Pair_AUC = DSN_info.Pair_AUC(1:100, 1:100);
    DSN_info.Gene_Name = DSN_info.Gene_Name(1:100);
end
if ~isequal(GE_Info.Gene_Name, DSN_info.Gene_Name), error(); end
Corr_mat = corr(zscore(GE_Info.Gene_Expression), 'Type', 'Spearman');
Corr_mat(1:n_gene+1:end) = 0;
Pair_Dist = (1-SetToZO(abs(Corr_mat))).^2;
clear Corr_mat
ind_auc = DSN_info.Pair_AUC(1:n_gene+1:end)';
pair_max = bsxfun(@max, ind_auc, ind_auc');
Pair_Dist = Pair_Dist + (1-SetToZO(DSN_info.Pair_AUC ./ pair_max)).^2;
ax_avg = bsxfun(@(x,y) (x+y)/2, ind_auc, ind_auc');
Pair_Dist = Pair_Dist + (1-SetToZO(ax_avg)).^2;
Net_Adj = -sqrt(Pair_Dist); % single(
clear DSN_info Pair_Dist ind_auc ax_avg pair_max
GE_Info.Gene_Expression = [];

% load fisheriris
% inds = ~strcmp(species,'setosa');
% Net_Adj = meas(inds,3:4);
% Node_Label = (strcmp(species(inds), 'virginica')==1)*2-1;
% Net_GeneName = GE_info.Gene_Name(1:n_gene);

% [~, sind] = sort(Net_Adj(:), 'Descend');
% [Pair_Info(:,1), Pair_Info(:,2)] = ind2sub([n_gene n_gene], sind(1:50));
% Pair_Info(Pair_Info(:,1)>Pair_Info(:,2), :) = [];
% a = GE_info.Gene_Name(Pair_Info);

%% Reduce Adj dimension
fprintf('Reducing dimension of Adj mat ...\n');
% [~, sid] = sort(std(Net_Adj), 'Descend');
[~, sid] = sort(abs(corr(Net_Adj, Census_isCancer, 'Type', 'Spearman')), 'Descend');
Net_Adj = Net_Adj(:, sid(1:2000));
% [coeff,score,latent,tsquared,explained,mu] = pca(Net_Adj, 'NumComponents', 10);
[coeff, ~, ~, ~, explained, ~] = pca(Net_Adj, 'NumComponents', 200);
Net_Adj = Net_Adj*coeff;

%% Run a cross-validation
Fold_Rank = zeros(n_gene, n_fold);
Model_AUC = zeros(n_fold, 1);
Fold_Index = crossvalind('KFold', Census_isCancer, n_fold);
for fi=1:n_fold
    %% Train SVM over adjacency matrix
    fprintf('Training SVM, Fold %d: ', fi);
    iTr = Fold_Index~=fi;
    SVMModel = SVM(Net_Adj(iTr, :), Census_isCancer(iTr), 'Linear', 0.5);
    if SVMModel.ClassNames(2)~=1, error(); end
    Model_AUC(fi,1) = getModelAUC(SVMModel, Net_Adj(~iTr, :), Census_isCancer(~iTr));
    fprintf('AUC is %0.2f\n', Model_AUC(fi)*100);
    
    %% Rerankong of genes
    [Pred_Lbl, Pred_prob] = predict(SVMModel, Net_Adj);
    %SVM_Confidence = max(Pred_prob, [], 2);
    %scatter(Net_Adj(:,1),Net_Adj(:,2), oscore(SVM_Confidence)*150+5, Node_Label*2);
    %colormap(lines(2));
    %scatter(Net_Adj(:,1),Net_Adj(:,2), oscore(Pred_prob(:,2))*150+5, Node_Label*2);
    Fold_Rank(:, fi) = Pred_prob(:,2);
end
NetWAS_Rank = SetToZO(mean(Fold_Rank, 2));
[~, NetWAS_Sind] = sort(NetWAS_Rank, 'Descend');

%% Comparing identified gene set
gwas_isHit = GWAS_Info{2} > -log10(0.005);
GWAS_Score = getAUC(Census_isCancer, gwas_isHit);

%% Plotting
n_TopList = floor(n_gene * linspace(0.001,0.1,50));
n_thresh = numel(n_TopList);
figure('Position', [100 100 1500 700]);
gwas_h = plot([0 n_thresh], GWAS_Score([1 1]), ':', 'LineWidth', 2, 'Color', [1 0 0]);
hold on
for ti=1:numel(n_TopList)
    Rank_arr = rand(n_gene, 1);
    Rank_arr(NetWAS_Sind(1:n_TopList)) = NetWAS_Rank(NetWAS_Sind(1:n_TopList));
    NetWAS_score = getAUC(Census_isCancer, Rank_arr);
    
    netwas_h(ti,1) = plot(ti, NetWAS_score, 'O', 'LineWidth', 1, 'Color', [0 0 1]);
end


%% Functions :::::::::::::::::::::::::::::::::::::
function SVMModel = SVM(X, Y, kernel_name, C, gamma)
if strcmpi(kernel_name, 'rbf')
    SVMModel = fitcsvm(X, Y, 'Verbose', 0, 'KernelFunction', kernel_name, 'BoxConstraint', C, 'KernelScale', gamma);
else
    SVMModel = fitcsvm(X, Y, 'Verbose', 1, 'KernelFunction', kernel_name, 'BoxConstraint', C); % 'Weights'
end
end

function [model_auc, pred_prob] = getModelAUC(SVMModel, xTe, lTe)
[~, pred_prob] = predict(SVMModel, xTe);
[~, pred_lbl] = max(pred_prob, [], 2);
model_auc = getAUC(lTe, pred_lbl, 50);
end

function data = SetToZO(data)
data = data - min(data(:));
data = data / max(data(:));
end