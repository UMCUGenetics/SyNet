
%% Initialization
clc;
clear;
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../_Utilities/');
addpath('../02_Perform_NOPs/');
n_lam = 20;

%% Load GE data
ge_name = getPath('SyNet');
fprintf('Loading expression from [%s] ...\n', ge_name);
load(ge_name, 'Gene_Expression', 'Patient_Label', 'Patient_Info', 'Study_Index');
Patient_Label = (Patient_Label==1)*2-1;

%% Get study index
si = 1;
iTr = Study_Index ~= si;
iTe = Study_Index == si;
if any(sum([iTr iTe],2)>1), error(); end

xTr = Gene_Expression(iTr,:);
xTe = Gene_Expression(iTe,:);
zTr = zscore(xTr);
zTe = zscore(xTe);
lTr = Patient_Label(iTr);
lTe = Patient_Label(iTe);

[opt_B, opt_fit] = lasso(xTr, lTr, 'NumLambda', n_lam, 'CV', 5);
% lasso_opt = {'lassoType', 't', 'CV', 5, 'relTol', 5e-1, 'n_lC', n_lam, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1};
% [opt_B, opt_fit] = lassoEx(zTr, lTr, lasso_opt{:});
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

plot([tr_auc_lam' te_auc_lam']);