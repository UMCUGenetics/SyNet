function S01_GenerateCV(Train_name, Train_Ind, Test_name, Test_Ind, rep_ind)
% S01_GenerateCV('ACES', 1:12, 'ACES', 1, 1)

%% Initialization
addpath('../_Utilities/');
[~,~] = mkdir('./CV_Files');
if strcmp(Train_name, Test_name) && any(ismember(Train_Ind, Test_Ind))
	fprintf('[i] Warning: Leakage found.\n');
	Train_Ind = setdiff(Train_Ind, Test_Ind);
	fprintf('Training set is cleaned to: %s\n', num2str(Train_Ind, '%2d, '));
	fprintf('Test set is: %s\n', num2str(Test_Ind, '%2d, '));
end

%% Load train information
tr_info.GEName = Train_name;
tr_info.GEPath = getPath(Train_name);
load(tr_info.GEPath, 'Patient_Label', 'Study_Index', 'Gene_Name');

n_sample = numel(Patient_Label);
TrnInd = find(ismember(Study_Index, Train_Ind));
n_Trn = numel(TrnInd);
sel_ind = randperm(n_Trn, floor(n_Trn*0.7));
tr_info.CVInd = false(n_sample, 1);
tr_info.CVInd(TrnInd(sel_ind)) = 1;
tr_info.iCvPar = Study_Index(tr_info.CVInd);
tr_info.Gene_Name = Gene_Name;
tr_info.Study_Ind = Train_Ind;

%% Load test information
te_info.GEName = Test_name;
te_info.GEPath = getPath(Test_name);
load(te_info.GEPath, 'Patient_Label', 'Study_Index', 'Gene_Name');
te_info.CVInd = ismember(Study_Index, Test_Ind);
te_info.Gene_Name = Gene_Name;
te_info.Study_Ind = Test_Ind;

%% Sanity check
if isequal(tr_info.GEPath, te_info.GEPath) && any(tr_info.CVInd & te_info.CVInd)
	error('Some samples test set are used in training.');
end

%% Saving
sav_name = sprintf('./CV_Files/CV_%s-%s_CVT02_Si%02d-Ri%03d.mat', Train_name, Test_name, Test_Ind, rep_ind);
if ~exist(sav_name, 'file')
	fprintf('Saving results in [%s]\n', sav_name);
	save(sav_name, 'tr_info', 'te_info', 'Patient_Label', 'Study_Index');
else
	error('[%s] CV file exists.\n', sav_name);
end
end
