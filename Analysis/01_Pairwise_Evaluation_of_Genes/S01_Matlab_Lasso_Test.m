%% Initialization
clc
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));

%% Load GE data
ge_name = getPath('SyNet');
fprintf('Loading expression from [%s] ...\n', ge_name);
load(ge_name, 'Gene_Expression', 'Patient_Label', 'Gene_Name', 'Study_Index');
Patient_Label = (Patient_Label==1)*2-1;

n_lam = 20;
si = 2;
iTr = Study_Index ~= si;
iTe = Study_Index == si;
xTr = Gene_Expression(iTr,:);
xTe = Gene_Expression(iTe,:);
zTr = zscore(xTr);
zTe = zscore(xTe);
lTr = Patient_Label(iTr);
lTe = Patient_Label(iTe);
n_zTr = sum(iTr);
n_zTe = sum(iTe);

[full_B, ~] = lasso(xTr, lTr, 'NumLambda', n_lam);
% lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1};
% [full_B, ~] = lassoEx(xTr, lTr, lasso_opt{:});

tr_pred = xTr*full_B;
te_pred = xTe*full_B;
auc_mat = zeros(2, n_lam);
for li=1:20
	[~, ~, ~, auc_mat(1,li)] = perfcurve(lTr, tr_pred(:,li), 1);
	[~, ~, ~, auc_mat(2,li)] = perfcurve(lTe, te_pred(:,li), 1);
end

close all
plot(auc_mat');