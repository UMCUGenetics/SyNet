function S01_Compare_Corr_in_test(Study_lst, Method_ind, Method_opt)
%{
%% Run
JOB_RANDOM_WAIT=1
PARAM=1000
for si in `seq 1 14`; do
	for mi in 3; do
		sbatch --job-name=Run-S${si}-M${mi}-$PARAM --mem=8GB --partition=general --qos=short --ntasks=1 --cpus-per-task=1 --time=4:00:00 --output=Logs/Log_Run_S${si}-M${mi}-$PARAM_%J_%a-%N.out run_Matlab.sh S01_Compare_Corr_in_test $si,$mi,$PARAM
	done
done
%}

%% Initialization
clc
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../01_Pairwise_Evaluation_of_Genes/');
addpath('../02_Perform_NOPs/');
if ismac
	fprintf('Warning: Running of debug mode ...\n');
	Study_lst = 1;
	Method_ind = 3;
end
if ~exist('Method_opt', 'var')
	Method_opt = 10000;
end

%% Load GE data
ge_name = getPath('SyNet');
fprintf('Loading expression from [%s] ...\n', ge_name);
load(ge_name, 'Gene_Expression', 'Patient_Label', 'Gene_Name', 'Study_Index');
% zData = zscore(Gene_Expression);
% [n_sample, n_gene] = size(Gene_Expression);
Patient_Label = (Patient_Label==1)*2-1;

%% Study loop
for si=1:numel(Study_lst)
	
	fprintf('Running for Study [%d]\n', Study_lst(si));
	%% Preparing data
	iTr = Study_Index ~= Study_lst(si);
	iTe = Study_Index == Study_lst(si);
	xTr = Gene_Expression(iTr,:);
	xTe = Gene_Expression(iTe,:);
	zTr = zscore(xTr);
	zTe = zscore(xTe);
	lTr = Patient_Label(iTr);
	lTe = Patient_Label(iTe);
	n_zTr = sum(iTr);
	n_zTe = sum(iTe);
	
	%% Method loop
	for mi=1:numel(Method_ind)
		fprintf('Running for Method [%d]\n', Method_ind(mi));
		
		%% Running Train procedure
		switch Method_ind(mi)
			case 1 % Lasso, 5CV
				lasso_opt = {'lassoType', 't', 'CV', 5, 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1};
				result = TrainLasso(zTr, lTr, zTe, lTe, lasso_opt);
				result.Method_Name = 'Single-5FCV';
			case 2
				lasso_opt = {'lassoType', 't', 'CV', 1, 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1, 'iCvPar', Study_Index(iTr)};
				result = TrainLasso(zTr, lTr, zTe, lTe, lasso_opt);
				result.Method_Name = 'Single-LOSO';
			case 3 % Lasso, Matlab implementation
				result = TrainMatlabLasso(zTr, lTr, zTe, lTe, Study_Index(iTr));
				result.Method_Name = 'Single-LOSOMat';
			case 4
				n_pair = Method_opt;
				Pair_Info = getTopPairs(Study_lst(si), n_pair);
				Gene_Set = unique(Pair_Info(:,1:2)', 'Stable');
				lasso_opt = {'lassoType', 't', 'CV', 1, 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1, 'iCvPar', Study_Index(iTr)};
				result = TrainLasso(zTr(:,Gene_Set), lTr, zTe(:,Gene_Set), lTe, lasso_opt);
				result.Gene_Set = Gene_Set;
				result.Method_Name = sprintf('SyNetGene-%05d-LOSO', n_pair);
			case 5
				n_pair = Method_opt;
				Pair_Info = getTopPairs(Study_lst(si), n_pair);
				mTr = zeros(n_zTr, n_pair);
				mTe = zeros(n_zTe, n_pair);
				[dTr, dTe] = CorrectGeneDirection(zTr, zTe, lTr);
				for pi=1:n_pair
					mTr(:,pi) = mean(dTr(:, Pair_Info(pi,1:2)), 2);
					mTe(:,pi) = mean(dTe(:, Pair_Info(pi,1:2)), 2);
				end
				
				lasso_opt = {'lassoType', 't', 'CV', 1, 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1, 'iCvPar', Study_Index(iTr)};
				result = TrainLasso(mTr, lTr, mTe, lTe, lasso_opt);
				result.Pair_Info = Pair_Info;
				result.Method_Name = sprintf('SyNetPair-%05d-LOSO', n_pair);
			otherwise
				error();
		end
		fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
		
		%% Saving results
		sav_name = sprintf('./Results/Res_Std-%02d_Met-%s.mat', Study_lst(si), result.Method_Name);
		fprintf('Saving results in [%s]\n', sav_name);
		save(sav_name, '-struct', 'result');
	end
end
end

%% ///////////// Functions
function result = TrainLasso(xTr, lTr, xTe, lTe, lasso_opt)
[opt_B, opt_fit] = lassoEx(xTr, lTr, lasso_opt{:});

vec_B = opt_B(:, opt_fit.IndexMinMSE);
g = xTr*vec_B;
tr_auc = getAUC(lTr, g, 50);

g = xTe*vec_B;
te_auc = getAUC(lTe, g, 50);

n_lam = size(opt_B, 2);
te_auc_lam = zeros(1, n_lam);
for i=1:n_lam
	te_auc_lam(i) = getAUC(lTe, xTe*opt_B(:,i), 50);
end

result.B = opt_B;
result.fit = opt_fit;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.te_auc_lam = te_auc_lam;
end

function result = TrainMatlabLasso(xTr, lTr, xTe, lTe, Fold_Index)
[n_sample, n_gene] = size(xTr);
n_lam = 20;
fprintf('Training over all samples [%d x %d]...\n', n_sample, n_gene);
[opt_B, opt_fit] = lasso(xTr, lTr, 'NumLambda', n_lam);
[~, ~, Fold_Index] = unique(Fold_Index, 'Stable');
n_iFold = max(Fold_Index);
fold_auc = zeros(n_iFold, n_lam);
fprintf('Training over folds ...\n');
for fi=1:n_iFold
	fprintf('Fold [%d/%d] ... \n', fi, n_iFold);
	[fold_B, ~] = lasso(xTr(Fold_Index~=fi,:), lTr(Fold_Index~=fi,:), 'NumLambda', n_lam);
	fold_pred = xTr(Fold_Index==fi,:)*fold_B;
	for li=1:20
		fold_auc(fi,li) = getAUC(lTr(Fold_Index==fi), fold_pred(:,li), 50);
	end
end
fold_mse = 1 - mean(fold_auc, 1);
[opt_fit.MinMSE, opt_fit.IndexMinMSE] = min(fold_mse);

vec_B = opt_B(:, opt_fit.IndexMinMSE);
g = xTr*vec_B;
tr_auc = getAUC(lTr, g, 50);

g = xTe*vec_B;
te_auc = getAUC(lTe, g, 50);

te_auc_lam = zeros(1, n_lam);
for i=1:n_lam
	te_auc_lam(i) = getAUC(lTe, xTe*opt_B(:,i), 50);
end

result.B = opt_B;
result.fit = opt_fit;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.te_auc_lam = te_auc_lam;
end

function Pair_Info = getTopPairs(Target_Study, n_pair)
%% Loading
dsn_path = '../01_Pairwise_Evaluation_of_Genes/Network_Files/';
dsn_name = sprintf('%sDSN_SyNetS%02d.mat', dsn_path, Target_Study);
fprintf('Loading DSN in [%s]\n', dsn_name);
dsn_data = load(dsn_name, 'Net_Adj');
n_gene = size(dsn_data.Net_Adj, 1);
if ~issymmetric(dsn_data.Net_Adj), error(); end

%% Sorting
[val, sid] = sort(dsn_data.Net_Adj(:), 'Descend');
sid = sid(1:n_pair*2);
val = val(1:n_pair*2);
Pair_Info = zeros(n_pair*2, 3);
[Pair_Info(:,1), Pair_Info(:,2)] = ind2sub([n_gene, n_gene], sid);
Pair_Info(:,3) = dsn_data.Net_Adj(sid);
if ~isequal(val, Pair_Info(:,3)), error(); end
Pair_Info(Pair_Info(:,1)>Pair_Info(:,2), :) = [];
end
%{
%% Prepare top pairs

%% Train using TP
mTr = zeros(n_zTr, n_pair);
mTe = zeros(n_zTe, n_pair);
[dTr, dTe] = CorrectGeneDirection(zTr, zTe, lTr);
for fi=1:n_pair
	mTr(:,fi) = mean(dTr(:, Pair_Info(fi,1:2)), 2);
	mTe(:,fi) = mean(dTe(:, Pair_Info(fi,1:2)), 2);
end
[tp_B, tp_fit] = lassoEx(mTr, lTr, lasso_opt{:});

tp_auc = zeros(1, 20);
for i=1:20
	tp_auc(i) = getAUC(lTe, mTe*tp_B(:,i), 50);
end

%% Plot
close all
figure('Position', [100 100 1400 500]);
hold on
plot(1:20, 1 - las_fit.MSE, 'o-', 'Color', [0.2 0.2 1.0]);
plot(las_fit.IndexMinMSE, las_auc(las_fit.IndexMinMSE), 'b*');
plot(1:20, las_auc, 'Color', [0 0 1]);

plot(1:20, 1 - tp_fit.MSE, 's-', 'Color', [1.0 0.2 0.2]);
plot(tp_fit.IndexMinMSE, tp_auc(tp_fit.IndexMinMSE), 'r*');
plot(1:20, tp_auc, 'Color', [1.0 0.0 0.0]);
legend('show');

%}
