function S02_Main_Code(Study_Ind, Repeat_Ind, Method_List, Method_Opt)
%{
for ri in `seq 4`; do
for si in `seq 1 12`; do
for np in 100 500 1000 5000 10000 20000; do
PARAM="$si,$ri,1:8,$np"
sbatch --job-name=JB-$PARAM --output=Logs/JB-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=6GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S02_Main_Code $PARAM
done
done
done
%}
clc;
if ismac
	fprintf('*** Warning!: Running on debug mode.\n');
	Method_List = 7;
	Method_Opt = 10000;
	Study_Ind = 1;
	Repeat_Ind = 1;
end

%% Initialization
SEED_INFO=rng; Run_ID=double(SEED_INFO.Seed);
fprintf('Run ID is: [%d]\n', Run_ID);
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath('../../../../Useful_Sample_Codes/getAUC');
addpath('../_Utilities/');
cv_path = './CV_Files/';
lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 0};

%% Load CV info
cv_id = sprintf('Si%02d_Ri%03d', Study_Ind, Repeat_Ind);
fprintf('[i] CV ID is defined to be [%s]\n', cv_id);
cv_info = dir([cv_path 'CV_*' cv_id '*.mat']);
cv_name = [cv_path cv_info.name];
fprintf('CV file found [%s]\n', cv_name);
cv_info = load(cv_name);

%% Load GE data
ge_name = getPath('SyNet');
fprintf('Loading expression from [%s] ...\n', ge_name);
GE_Data = load(ge_name, 'Gene_Expression', 'Patient_Label', 'Patient_Info', 'Study_Index', 'Study_Name', 'Gene_Name');
GE_Data.Patient_Label = (GE_Data.Patient_Label==1)*2-1;

%% Prepare data
iTr = cv_info.cv_obj.iTr;
iTe = cv_info.cv_obj.iTe;
if any(sum([iTr iTe], 2)>1), error(); end
Dataset_info.zTr = zscore(GE_Data.Gene_Expression(iTr,:));
Dataset_info.zTe = zscore(GE_Data.Gene_Expression(iTe,:));
Dataset_info.lTr = GE_Data.Patient_Label(iTr);
Dataset_info.lTe = GE_Data.Patient_Label(iTe);
Dataset_info.Gene_Name = GE_Data.Gene_Name;
[~, ~, Dataset_info.Fold_Index] = unique(GE_Data.Study_Index(iTr), 'Stable');

%% Main Loop
initialVars = [who; 'mi'; 'initialVars'];
for mi=1:numel(Method_List)
	clearvars('-except', initialVars{:});
	fprintf('Running method [%d]...\n', Method_List(mi));
	switch Method_List(mi)
		case 1
			Method_Name = 'NMC';
			fprintf('/// Method name is [%s]\n', Method_Name);
			result = Perf_NMC(Dataset_info);
		case 2
			Method_Name = 'LEx';
			fprintf('/// Method name is [%s]\n', Method_Name);
			result = LassoWithCV(@lassoEx, Dataset_info.zTr, Dataset_info.lTr, Dataset_info.zTe, Dataset_info.lTe, Dataset_info.Fold_Index, lasso_opt);
		case 3
			n_pair = Method_Opt;
			Method_Name = sprintf('NMC-DSN-%05d', n_pair);
			fprintf('/// Method name is [%s]\n', Method_Name);
			CMB_Data = Dataset_info;
			
			Pair_Info = getTopPairs(cv_info.cv_obj.Test_Ind, n_pair, GE_Data.Study_Index);
			Gene_Set = unique(Pair_Info(:,1:2)', 'Stable');
			CMB_Data.zTr = CMB_Data.zTr(:, Gene_Set);
			CMB_Data.zTe = CMB_Data.zTe(:, Gene_Set);
			
			result = Perf_NMC(CMB_Data);
			result.n_pair = n_pair;
		case 4
			n_pair = Method_Opt;
			Method_Name = sprintf('LEx-DSN-%05d', n_pair);
			fprintf('/// Method name is [%s]\n', Method_Name);
			CMB_Data = Dataset_info;
			
			Pair_Info = getTopPairs(cv_info.cv_obj.Test_Ind, n_pair, GE_Data.Study_Index);
			Gene_Set = unique(Pair_Info(:,1:2)', 'Stable');
			CMB_Data.zTr = CMB_Data.zTr(:, Gene_Set);
			CMB_Data.zTe = CMB_Data.zTe(:, Gene_Set);
			
			result = LassoWithCV(@lassoEx, CMB_Data.zTr, CMB_Data.lTr, CMB_Data.zTe, CMB_Data.lTe, CMB_Data.Fold_Index, lasso_opt);
			result.n_pair = n_pair;
		case 5
			n_pair = Method_Opt;
			Method_Name = sprintf('NMC-DsnGrp-%05d', n_pair);
			fprintf('/// Method name is [%s]\n', Method_Name);
			CMB_Data = Dataset_info;
			
			Pair_Info = getTopPairs(cv_info.cv_obj.Test_Ind, n_pair, GE_Data.Study_Index);
			
			Nei_grp = arrayfun(@(i) Pair_Info(i,1:2), (1:n_pair)', 'UniformOutput', false);
			[CMB_Data.zTr, CMB_Data.zTe] = CombineGenes(Nei_grp, CMB_Data.zTr, CMB_Data.zTe, 2);
			
			result = Perf_NMC(CMB_Data);
			result.n_pair = n_pair;
		case 6
			n_pair = Method_Opt;
			Method_Name = sprintf('LEx-DsnGrp-%05d', n_pair);
			fprintf('/// Method name is [%s]\n', Method_Name);
			CMB_Data = Dataset_info;
			
			Pair_Info = getTopPairs(cv_info.cv_obj.Test_Ind, n_pair, GE_Data.Study_Index);
			
			Nei_grp = arrayfun(@(i) Pair_Info(i,1:2), (1:n_pair)', 'UniformOutput', false);
			[CMB_Data.zTr, CMB_Data.zTe] = CombineGenes(Nei_grp, CMB_Data.zTr, CMB_Data.zTe, 2);
			
			result = LassoWithCV(@lassoEx, CMB_Data.zTr, CMB_Data.lTr, CMB_Data.zTe, CMB_Data.lTe, CMB_Data.Fold_Index, lasso_opt);
			result.n_pair = n_pair;
		case 7
			n_pair = Method_Opt;
			Method_Name = sprintf('LEx-DsnIReg-%05d', n_pair);
			fprintf('/// Method name is [%s]\n', Method_Name);
			CMB_Data = Dataset_info;
			
			Pair_Info = getTopPairs(cv_info.cv_obj.Test_Ind, n_pair, GE_Data.Study_Index);
			
			Nei_grp = arrayfun(@(i) Pair_Info(i,1:2), (1:n_pair)', 'UniformOutput', false);
			[CMB_Data.zTr, CMB_Data.zTe] = IntegGenes(Nei_grp, CMB_Data.zTr, CMB_Data.lTr, CMB_Data.zTe, CMB_Data.Fold_Index);
			
			result = LassoWithCV(@lassoEx, CMB_Data.zTr, CMB_Data.lTr, CMB_Data.zTe, CMB_Data.lTe, CMB_Data.Fold_Index, lasso_opt);
			result.n_pair = n_pair;
		case 8
			Method_Name = sprintf('LEx-FoldReg');
			fprintf('/// Method name is [%s]\n', Method_Name);
			
			result = LassoExOverFolds(Dataset_info.zTr, Dataset_info.lTr, Dataset_info.zTe, Dataset_info.lTe, Dataset_info.Fold_Index);
	end
	
	%% Storing in variable
	result.Method_Name = Method_Name;
	result.Study_Ind = Study_Ind;
	result.Repeat_Ind = Repeat_Ind;
	result.CV_Info = cv_info;
	result.CV_ID = cv_id;
	fprintf('@@@@@ Final test performance for [%s] is [%0.2f%%] AUC.\n', result.Method_Name, result.te_auc*100);
	
	%% Saving results
	sav_path = './Results/';
	[~, ~] = mkdir(sav_path);
	sav_name = sprintf('%sRes_%s_%s_%07d.mat', sav_path, result.Method_Name, result.CV_ID, Run_ID);
	fprintf('Saving in [%s]\n', sav_name);
	save(sav_name, '-struct', 'result');
	
	fprintf('*************\n\n');
end

