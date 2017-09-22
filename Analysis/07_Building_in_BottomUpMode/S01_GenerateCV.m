function S01_GenerateCV(data_name, Train_Ind, Test_Ind, rep_ind)
% S01_GenerateCV('SyNet', 1:14, 1, 101)

%% Initialization
addpath('../_Utilities/');
[~,~] = mkdir('./CV_Files');
if any(ismember(Train_Ind, Test_Ind))
	fprintf('[i] Warning: Leakage found. Training set is cleaned to: ');
	Train_Ind = setdiff(Train_Ind, Test_Ind);
	fprintf('%d,', Train_Ind);
	fprintf('\n');
end

%% Load Study information
GeneExpression_Path = getPath(data_name);
load(GeneExpression_Path, 'Patient_Label', 'Study_Index');
n_sample = numel(Patient_Label);

%% Main code
cv_obj.iTe = ismember(Study_Index, Test_Ind);

TrnInd = find(ismember(Study_Index, Train_Ind));
n_Trn = numel(TrnInd);
sel_ind = randperm(n_Trn, floor(n_Trn*0.7));
cv_obj.iTr = false(n_sample, 1);
cv_obj.iTr(TrnInd(sel_ind)) = 1;
if any(sum([cv_obj.iTr cv_obj.iTe],2)>1), error(); end

cv_obj.Train_Ind = Train_Ind;
cv_obj.Test_Ind = Test_Ind;

%% Saving
sav_name = sprintf('./CV_Files/CV_%s_CVT01_Si%02d_Ri%03d.mat', data_name, Test_Ind, rep_ind);
if ~exist(sav_name, 'file')
	fprintf('Saving results in [%s]\n', sav_name);
	save(sav_name, 'cv_obj', 'Patient_Label', 'Study_Index');
else
	error('[%s] CV file exists.\n', sav_name);
end
end
