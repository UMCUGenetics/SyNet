function S04_Predict_StudyLabel(Train_si, Test_si)
%% Run
%{
for i in `seq 1 12`; do
PARAM=13,$i; sbatch --job-name=Study-${PARAM} --mem=6GB --partition=general --qos=short --ntasks=1 --cpus-per-task=1 --time=4:00:00 --output=Logs/Log_Study_${PARAM}_%J_%a-%N.out run_Matlab.sh S04_Predict_StudyLabel ${PARAM}
done
%}

%% Initialization
clc % S04_Predict_StudyLabel(1:12,11)
SEED_INFO=rng; Run_ID=double(SEED_INFO.Seed);
fprintf('Run ID is: [%d]\n', Run_ID);
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../_Utilities/');
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
GE_Data.Patient_Label = (GE_Data.Study_Index==Test_si)*2-1;
n_sample = size(GE_Data.Gene_Expression, 1);
cv_obj = cvpartition(GE_Data.Patient_Label, 'kFold', 5);
GE_Data.Fold_Index = zeros(n_sample, 1);
for fi=1:cv_obj.NumTestSets
	GE_Data.Fold_Index(cv_obj.test(fi), 1) = fi;
end

%% Set names
if numel(Train_si)>1
	in_Trn = ismember(GE_Data.Study_Index, Train_si);
	Train_Set = unique(GE_Data.Patient_Info.Source_Study(in_Trn));
	TrName = strjoin(Train_Set, '-');
else
	TrName = GE_Data.Study_Name{Train_si};
end
TeName = strjoin(GE_Data.Study_Name(Test_si),',');

%% Prepare the sample indices
fprintf('Preparing data for [%s, %s]\n', TrName, TeName);
iTr = ismember(GE_Data.Fold_Index, 1:4);
iTe = ismember(GE_Data.Fold_Index, 5);
if any(sum([iTr iTe],2)>1), error(); end
fprintf('#Tr=%d, #Te=%d\n', sum(iTr), sum(iTe));

%% Get study index
xTr = GE_Data.Gene_Expression(iTr,:);
xTe = GE_Data.Gene_Expression(iTe,:);
zTr = zscore(xTr);
zTe = zscore(xTe);
lTr = GE_Data.Patient_Label(iTr);
lTe = GE_Data.Patient_Label(iTe);

%% Training the classifier
Classifier_Name = 'LassoEx';
func = @lassoEx;
lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 0};
fprintf('Running [%s] ...\n', Classifier_Name);
result = LassoWithCV(func, zTr, lTr, zTe, lTe, GE_Data.Fold_Index(iTr), lasso_opt);

%% Storing in variable
result.TrName = TrName;
result.TeName = TeName;
result.Train_si = Train_si;
result.Test_si = Test_si;
result.Study_lst = GE_Data.Study_Name;
result.Classifier_Name = Classifier_Name;
result.iTr = iTr;
result.iTe = iTe;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);

%% Saving results
sav_path = './S04_Predict_StudyLabel/';
[~, ~] = mkdir(sav_path);
sav_name = sprintf('%sRes_%s_%s_%s_%07d.mat', sav_path, Classifier_Name, result.TrName, result.TeName, Run_ID);
fprintf('Saving in [%s]\n', sav_name);
save(sav_name, '-struct', 'result');
end

