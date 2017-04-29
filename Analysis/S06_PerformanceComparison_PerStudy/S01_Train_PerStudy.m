function S01_Train_PerStudy(Train_si, Test_si)
%% Run
%{
for PARAM in 10:14,11 8:13,13; do
sbatch --job-name=Run-${PARAM} --mem=8GB --partition=general --qos=short --ntasks=1 --cpus-per-task=1 --time=4:00:00 --output=Logs/Log_Run_${PARAM}_%J_%a-%N.out run_Matlab.sh S01_Train_PerStudy ${PARAM}
done
%}

%% Initialization
clc % S01_Train_PerStudy([1:12 14],13)
SEED_INFO=rng; Run_ID=double(SEED_INFO.Seed);
fprintf('Run ID is: [%d]\n', Run_ID);
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../_Utilities/');
Ratio_lst = [0.7]; %[0.3 0.5 0.7 0.9 1.0];
if any(ismember(Train_si, Test_si))
	fprintf('[i] Warning: Leakage found. Training set is cleaned to: ');
	Train_si = setdiff(Train_si, Test_si);
	fprintf('%d,', Train_si);
	fprintf('\n');
end

%% Load GE data
ge_name = getPath('SyNet');
fprintf('Loading expression from [%s] ...\n', ge_name);
GE_Data = load(ge_name, 'Gene_Expression', 'Patient_Label', 'Patient_Info', 'Study_Index', 'Study_Name');
GE_Data.Patient_Label = (GE_Data.Patient_Label==1)*2-1;
n_sample = size(GE_Data.Gene_Expression, 1);

%% Main Loop
for ri=1:numel(Ratio_lst)
	%% Choose ratio
	Sample_Ratio = Ratio_lst(ri);
	fprintf('Ratio is set to [%0.2f]\n', Sample_Ratio);
	
	%% Prepare the sample indices
	if numel(Train_si)>1
		in_Trn = ismember(GE_Data.Study_Index, Train_si);
		Train_Set = unique(GE_Data.Patient_Info.Source_Study(in_Trn));
		TrName = strjoin(Train_Set, '-');
	else
		TrName = GE_Data.Study_Name{Train_si};
	end
	TeName = strjoin(GE_Data.Study_Name(Test_si),',');
	fprintf('Preparing data for [%s, %s]\n', TrName, TeName);
	TrnInd = find(ismember(GE_Data.Study_Index, Train_si));
	n_iTr = numel(TrnInd);
	sel_ind = randperm(n_iTr, floor(n_iTr*Sample_Ratio));
	iTr = false(n_sample, 1);
	iTr(TrnInd(sel_ind)) = 1;
	iTe = ismember(GE_Data.Study_Index, Test_si);
	if any(sum([iTr iTe],2)>1), error(); end
	fprintf('#Tr=%d, #Te=%d\n', sum(iTr), sum(iTe));
	
	if ~exist('initialVars', 'var')
		initialVars = [who; 'ci'; 'initialVars'];
	end
	for ci=1:1
		clearvars('-except', initialVars{:});
		
		%% Choose classifier
		switch ci
			case 1
				Classifier_Name = 'LassoEx';
				func = @lassoEx;
				lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 0};
			case 2
				Classifier_Name = 'Lasso';
				func = @lasso;
				lasso_opt = {'NumLambda', 20};
		end
		
		%% Get study index
		xTr = GE_Data.Gene_Expression(iTr,:);
		xTe = GE_Data.Gene_Expression(iTe,:);
		zTr = zscore(xTr);
		zTe = zscore(xTe);
		lTr = GE_Data.Patient_Label(iTr);
		lTe = GE_Data.Patient_Label(iTe);
		
		%% Training the classifier
		fprintf('Running [%s] ...\n', Classifier_Name);
		result = LassoWithCV(func, zTr, lTr, zTe, lTe, GE_Data.Study_Index(iTr), lasso_opt);
		
		%% Storing in variable
		result.TrName = TrName;
		result.TeName = TeName;
		result.Train_si = Train_si;
		result.Test_si = Test_si;
		result.Study_lst = GE_Data.Study_Name;
		result.Classifier_Name = Classifier_Name;
		result.Sample_Ratio = Sample_Ratio;
		result.iTr = iTr;
		result.iTe = iTe;
		fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);
		
		%% Saving results
		sav_path = './S01_Performance_PerStudy/';
		[~, ~] = mkdir(sav_path);
		sav_name = sprintf('%sRes_%s_%s_%s_SR%0.2f_%07d.mat', sav_path, Classifier_Name, result.TrName, result.TeName, Sample_Ratio, Run_ID);
		fprintf('Saving in [%s]\n', sav_name);
		save(sav_name, '-struct', 'result');
	end
end
end

