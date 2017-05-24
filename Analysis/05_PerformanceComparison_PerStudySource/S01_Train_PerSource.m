function S01_Train_PerSource(Train_si, Test_si, Classifier_Ind)
%% Run
%{
for i in 1,2 1,3 2,3 2,1 3,1 3,2; do
PARAM=$i,1; sbatch --job-name=Run-${PARAM} --mem=7GB --partition=general --qos=short --ntasks=1 --cpus-per-task=1 --time=4:00:00 --output=Logs/Log_Run_${PARAM}_%J_%a-%N.out run_Matlab.sh S01_Train_PerSource ${PARAM}
done
%}

%% Initialization
clc % S01_Train_PerSource(2,1,1)
SEED_INFO=rng; Run_ID=double(SEED_INFO.Seed);
fprintf('Run ID is: [%d]\n', Run_ID);
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../_Utilities/');
n_lam = 20;
Ratio_lst = [0.3 0.5 0.7 0.9 1.0];
switch Classifier_Ind
	case 1
		Classifier_Name = 'Lasso';
	case 2
		Classifier_Name = 'LassoEx';
end

%% Load GE data
ge_name = getPath('SyNet');
fprintf('Loading expression from [%s] ...\n', ge_name);
load(ge_name, 'Gene_Expression', 'Patient_Label', 'Patient_Info');
Patient_Label = (Patient_Label==1)*2-1;
[n_sample, n_gene] = size(Gene_Expression);

%% Main Loop
initialVars = [who; 'ri'; 'initialVars'];
for ri=1:numel(Ratio_lst)
	clearvars('-except', initialVars{:});
	
	Sample_Ratio = Ratio_lst(ri);
	fprintf('Ratio is set to [%0.2f]\n', Sample_Ratio);
	
	%% Get study index
	[Study_lst, ~, Study_Index] = unique(Patient_Info.Source_Study);
	fprintf('Preparing data for [%s, %s]\n', Study_lst{Train_si}, Study_lst{Test_si});
	result.TrName = Study_lst{Train_si};
	result.TeName = Study_lst{Test_si};
	
	iTe = Study_Index == Test_si;
	TrnInd = find(Study_Index == Train_si);
	n_iTr = numel(TrnInd);
	sel_ind = randperm(n_iTr, floor(n_iTr*Sample_Ratio));
	iTr = false(n_sample, 1);
	iTr(TrnInd(sel_ind)) = 1;
	if any(sum([iTr iTe],2)>1), error(); end
	fprintf('#Tr=%d, #Te=%d\n', sum(iTr), sum(iTe));
	
	xTr = Gene_Expression(iTr,:);
	xTe = Gene_Expression(iTe,:);
	zTr = zscore(xTr);
	zTe = zscore(xTe);
	lTr = Patient_Label(iTr);
	lTe = Patient_Label(iTe);
	
	%% Training the classifier
	fprintf('Running [%s] ...\n', Classifier_Name);
	switch Classifier_Name
		case 'Lasso'
			[opt_B, opt_fit] = lasso(zTr, lTr, 'NumLambda', n_lam, 'CV', 5);
		case 'LassoEx'
			lasso_opt = {'lassoType', 't', 'CV', 5, 'relTol', 5e-1, 'n_lC', n_lam, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1};
			[opt_B, opt_fit] = lassoEx(zTr, lTr, lasso_opt{:});
	end
	n_lam = size(opt_B, 2);
	
	%% Evaluation
	vec_B = opt_B(:, opt_fit.IndexMinMSE);
	g = zTr*vec_B;
	tr_auc = getAUC(lTr, g, 50);
	
	g = zTe*vec_B;
	te_auc = getAUC(lTe, g, 50);
	
	tr_auc_lam = zeros(1, n_lam);
	te_auc_lam = zeros(1, n_lam);
	for i=1:n_lam
		tr_auc_lam(i) = getAUC(lTr, zTr*opt_B(:,i), 50);
		te_auc_lam(i) = getAUC(lTe, zTe*opt_B(:,i), 50);
	end
	
	%% Storing in variable
	result.Train_si = Train_si;
	result.Test_si = Test_si;
	result.Study_lst = Study_lst;
	result.Classifier_Name = Classifier_Name;
	result.Sample_Ratio = Sample_Ratio;
	result.iTr = iTr;
	result.iTe = iTe;
	result.B = opt_B;
	result.fit = opt_fit;
	result.tr_auc = tr_auc;
	result.te_auc = te_auc;
	result.tr_auc_lam = tr_auc_lam;
	result.te_auc_lam = te_auc_lam;
	fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
	
	%% Saving results
	sav_path = './S01_Performance_PerSource/';
	[~, ~] = mkdir(sav_path);
	sav_name = sprintf('%sRes_%s_%s_%s_SR%0.2f_%07d.mat', sav_path, Classifier_Name, Study_lst{Train_si}, Study_lst{Test_si}, Sample_Ratio, Run_ID);
	fprintf('Saving in [%s]\n', sav_name);
	save(sav_name, '-struct', 'result');
end
end

