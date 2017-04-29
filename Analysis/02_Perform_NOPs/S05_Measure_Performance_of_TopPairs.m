clc;
clear;

%% Initialization
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
dsn_path = '../01_Pairwise_Evaluation_of_Genes/Network_Files/';
study_lst = 1:14;
n_study = numel(study_lst);
tresh_lst = [50 100 200 500 1000 2000 5000 10000 20000];
n_tresh = numel(tresh_lst);
n_pair = max(tresh_lst);
poolobj = gcp('nocreate'); % If no pool, do not create new one.
if ~ismac && isempty(poolobj)
	fprintf('Opening workers...\n');
    parpool(15);
end

%% Load GE data
GE_data = load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Expression', 'Patient_Label', 'Gene_Name', 'Study_Index');
GE_data.Patient_Label = (GE_data.Patient_Label==1)*2-1;

%% Main loop
auc_lst = cell(n_study, 1);
for si=1:n_study
	%% Load expression
	iTr = GE_data.Study_Index ~= study_lst(si);
	iTe = GE_data.Study_Index == study_lst(si);
	lTr = GE_data.Patient_Label(iTr);
	lTe = GE_data.Patient_Label(iTe);
	xTr = GE_data.Gene_Expression(iTr,:);
	xTe = GE_data.Gene_Expression(iTe,:);
	
	%% Load top pairs
	dsn_name = sprintf('%sDSN_SyNetS%02d.mat', dsn_path, study_lst(si));
	fprintf('Loading DSN in [%s]\n', dsn_name);
	dsn_data = load(dsn_name, 'Net_Adj', 'Gene_Name'); 
	if ~isequal(dsn_data.Gene_Name, GE_data.Gene_Name), error(); end
	n_gene = numel(dsn_data.Gene_Name);
	[val, sid] = sort(dsn_data.Net_Adj(:), 'Descend');
	sid = sid(1:n_pair);
	val = val(1:n_pair);
	Pair_Info = zeros(n_pair, 3);
	[Pair_Info(:,1), Pair_Info(:,2)] = ind2sub([n_gene, n_gene], sid);
	Pair_Info(:,3) = dsn_data.Net_Adj(sid);
	if ~isequal(val, Pair_Info(:,3)), error(); end
	dsn_data = [];
	
	%% Loop over top pairs
	indv_B = [];
	for ti=1:n_tresh
		%% Select pairs
		g_set = unique(Pair_Info(1:tresh_lst(ti),1:2), 'Stable');
		auc_lst{si}(ti,1) = TrainLasso(xTr(:,g_set), lTr, xTe(:,g_set), lTe);
		
		%% Train classical lasso
		if isempty(indv_B)
			[~, ~, indv_B] = TrainLasso(xTr, lTr, xTe, lTe);
		end
		[val, sid] = sort(abs(indv_B), 'Descend');
		vec_B = indv_B;
		vec_B(sid(numel(g_set)+1:end)) = 0;
		auc_lst{si}(ti,2) = getAUC(lTe, xTe*vec_B, 50);
		
		%% Select random genes
		r_set = randperm(n_gene, numel(g_set));
		auc_lst{si}(ti,3) = TrainLasso(xTr(:,r_set), lTr, xTe(:,r_set), lTe);
	end
	
	%% Plotting
	% plot(auc_lst);
end
method_lst = {'TopPairs', 'All', 'Random'};
save('./S05_Individual_vs_TopPairs_Performance.mat', 'auc_lst', 'tresh_lst', 'method_lst');

function [te_auc, tr_auc, vec_B] = TrainLasso(xTr, lTr, xTe, lTe)

%% Initialization
lasso_opt = {'lassoType', 't', 'CV', 5, 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 0};
iCvPar = cvpartition(lTr, 'Kfold', 5);
n_feat = size(xTr, 2);
zTr = zscore(xTr);
zTe = zscore(xTe);

%% Traning the final model
fprintf('Training the final model over [%d] features...\n', n_feat);
[opt_B, opt_fit] = lassoEx(zTr, lTr, lasso_opt{:}, 'iCvPar', iCvPar);
fprintf('Final training is done. [%d] non-zero features identified.\n', sum(abs(opt_B(:, opt_fit.IndexMinMSE))>0));

%% Evaluating the model
vec_B = opt_B(:, opt_fit.IndexMinMSE);
% [~, topB_ind] = sort(abs(vec_B), 'Descend');
% vec_B(topB_ind(MAX_N_SUBNET+1:end)) = 0;

%% Plotting
% lassoPlot(opt_B, opt_fit, 'PlotType', 'CV')
% figure();
% plot(fliplr(sum(opt_B>0)));

g = zTr*vec_B;
tr_auc = getAUC(lTr, g, 50);

g = zTe*vec_B;
te_auc = getAUC(lTe, g, 50);
end