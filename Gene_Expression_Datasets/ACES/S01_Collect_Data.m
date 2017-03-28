%% Download data
% ! wget http://ccb.nki.nl/software/aces/ACES.tar.gz
% ! tar -xvzf ACES.tar.gz

%% Initialization
clc;
clear;
h5f_All = './ACES/experiments/data/U133A_combat.h5';
h5f_PerStudy = './ACES/experiments/data/U133A_combat_separate_studies.h5';

%% Loading full data
fprintf('Reading full dataset.\n');
Expression_Data = h5read(h5f_All, '/U133A_combat_RFS/ExpressionData')';
Gene_Entrez = h5read(h5f_All, '/U133A_combat_RFS/GeneLabels');
Patient_Class_Labels = h5read(h5f_All, '/U133A_combat_RFS/PatientClassLabels');
Patient_ID = h5read(h5f_All, '/U133A_combat_RFS/PatientLabels');
[n_pat, n_gene] = size(Expression_Data);
Pat_Map = containers.Map(Patient_ID, 1:n_pat);
if n_pat~=Pat_Map.Count, error(); end

%% Loading per-study data
step = 1;
Study_Name = {};
Study_Index = zeros(n_pat, 1);
h_info = h5info(h5f_PerStudy);
for gi=1:numel(h_info.Groups)
	grp_name = h_info.Groups(gi).Name;
	if ~contains(grp_name, '_RFS_')
		continue;
	end
	fprintf('Reading [%s] group.\n', grp_name);
	study_info = strsplit(grp_name, '_');
	Study_Name{step, 1} = study_info{end};
	
	Batch_Expr = h5read(h5f_PerStudy, [grp_name '/ExpressionData'])';
	Batch_Entz = h5read(h5f_PerStudy, [grp_name '/GeneLabels']);
	Batch_PatLbl = h5read(h5f_PerStudy, [grp_name '/PatientClassLabels']);
	Batch_PatID = h5read(h5f_PerStudy, [grp_name '/PatientLabels']);
	if ~isequal(Batch_Entz, Gene_Entrez), error(); end
	
	%% Check Expressions
	for pi=1:numel(Batch_PatID)
		pat_ind = Pat_Map(Batch_PatID{pi});
		if ~isequal(Batch_Expr(pi,:), Expression_Data(pat_ind,:)) || ~isequal(Batch_PatID{pi}, Patient_ID{pat_ind})
			error();
		end
		Study_Index(pat_ind) = step;
	end
	step = step + 1;
end

%% Run on python
% import cPickle
% from scipy.io import savemat
% fid = open('./ACES/experiments/data/U133A_combat_classlabels.pickle')
% ClinicalFeatures = cPickle.load(fid)
% fid.close()
% savemat('./ACES_ClinicalFeatures.mat', {'ClinicalFeatures':ClinicalFeatures})

%% Load clinical variables
load('./ACES_ClinicalFeatures.mat');
Clic_StudyGSM = deblank(cellstr(ClinicalFeatures{1}));
Clic_PatID = deblank(cellstr(ClinicalFeatures{2}));
Clic_tmp = deblank(cellstr(ClinicalFeatures{3}));
Clic_CancerType = deblank(cellstr(ClinicalFeatures{8,1}{1}));
for pi=1:size(PatientLabels,1)
    index = find(ismember(Clic_PatID, deblank(PatientLabels{pi})));
    if numel(index)>1
        warning('Duplicate samples detected: ');
        fprintf('%d, ', index);
        fprintf('\n');
    elseif numel(index)==0
        error('Could not find the ID.\n')
    else
        CancerTypeName{pi,1} = Clic_CancerType{index};
    end
end
[CancerTypeList, ~, CancerTypeIndex]=unique(CancerTypeName);
