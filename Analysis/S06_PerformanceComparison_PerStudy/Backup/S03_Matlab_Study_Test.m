function S03_Matlab_Study_Test(Train_si, Test_si)
%% Run
% si=1;sj=2;sbatch --job-name=Test-S${si}-S${sj} --mem=8GB --partition=general --qos=short --ntasks=1 --cpus-per-task=1 --time=4:00:00 --output=Logs/Log_Test_S${si}_S${sj}_%J_%a-%N.out run_Matlab.sh S03_Matlab_Study_Test $si,$sj

%% Initialization
clc % Train_si = 1; Test_si = 2;
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../01_Pairwise_Evaluation_of_Genes/');
addpath('../02_Perform_NOPs/');
Classifier_Name='Lasso';
% Classifier_Name='LassoEx';
n_lam = 20;

%% Load GE data
ge_name = getPath('SyNet');
fprintf('Loading expression from [%s] ...\n', ge_name);
load(ge_name, 'Gene_Expression', 'Patient_Label', 'Patient_Info');
Patient_Label = (Patient_Label==1)*2-1;

%% Get study index
[Study_lst, ~, Study_Index] = unique(Patient_Info.Source_Study);
fprintf('Preparing data for [%s, %s]\n', Study_lst{Train_si}, Study_lst{Test_si});
result.TrName = Study_lst{Train_si};
result.TeName = Study_lst{Test_si};
iTr = Study_Index == Train_si;
iTe = Study_Index == Test_si;
if any(sum([iTr iTe],2)>1), error(); end

xTr = Gene_Expression(iTr,:);
xTe = Gene_Expression(iTe,:);
zTr = zscore(xTr);
zTe = zscore(xTe);
lTr = Patient_Label(iTr);
lTe = Patient_Label(iTe);

fprintf('Running [%s] ...\n', Classifier_Name);
if strcmp(Classifier_Name, 'Lasso')
	[opt_B, opt_fit] = lasso(xTr, lTr, 'NumLambda', n_lam, 'CV', 5);
else
	lasso_opt = {'lassoType', 't', 'CV', 5, 'relTol', 5e-1, 'n_lC', n_lam, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1};
	[opt_B, opt_fit] = lassoEx(zTr, lTr, lasso_opt{:});
end
n_lam = size(opt_B, 2);

vec_B = opt_B(:, opt_fit.IndexMinMSE);
g = zTr*vec_B;
tr_auc = getAUC(lTr, g, 50);

g = zTe*vec_B;
te_auc = getAUC(lTe, g, 50);

tr_auc_lam = zeros(1, n_lam);
te_auc_lam = zeros(1, n_lam);
for i=1:n_lam
	tr_auc_lam(i) = getAUC(lTr, xTr*opt_B(:,i), 50);
	te_auc_lam(i) = getAUC(lTe, xTe*opt_B(:,i), 50);
end

result.B = opt_B;
result.fit = opt_fit;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
result.tr_auc_lam = tr_auc_lam;
result.te_auc_lam = te_auc_lam;
fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', result.te_auc*100);

%% Saving results
sav_path = './S01_Matlab_Study_Test';
[~, ~] = mkdir(sav_path);
sav_name = sprintf('%s/%s_%s_%s.mat', sav_path, Classifier_Name, Study_lst{Train_si}, Study_lst{Test_si});
fprintf('Saving in [%s]\n', sav_name);
save(sav_name, '-struct', 'result');
end

