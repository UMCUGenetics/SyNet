function S08_GenerateCV_ForOneStudy(GE_name, Target_Name, rep_ind)
% S08_GenerateCV_ForOneStudy('SyNet', 'TCGA', 1)

%% Initialization
addpath('../_Utilities/');

%% Load Study data
tr_info.GEName = GE_name;
tr_info.GEPath = getPath(GE_name);
load(tr_info.GEPath, 'Patient_Label', 'Study_Index', 'Gene_Name', 'Study_Name');
n_sample = numel(Patient_Label);
tr_info.CVInd = false(n_sample, 1);
tr_info.Gene_Name = Gene_Name;
te_info = tr_info;

%% Make CV indices
Study_Ind = find(strcmp(Study_Name, Target_Name));
is_val = Study_Index==Study_Ind;
Fold_Index = crossvalind('KFold', Patient_Label(is_val), 10);
tr_info.CVInd(is_val) = Fold_Index  > 3;
te_info.CVInd(is_val) = Fold_Index <= 3;

%% Sanity check
if isequal(tr_info.GEPath, te_info.GEPath) && any(tr_info.CVInd & te_info.CVInd)
	error('Some samples test set are used in training.');
end

%% Saving
sav_name = sprintf('./CV_Files/CV_%s-%s_CVT51_Si%02d-Ri%03d.mat', GE_name, GE_name, Study_Ind, rep_ind);
if ~exist(sav_name, 'file')
	fprintf('Saving results in [%s]\n', sav_name);
	save(sav_name, 'tr_info', 'te_info', 'Patient_Label', 'Study_Index');
else
	error('[%s] CV file exists.\n', sav_name);
end
end
