%% Download data
% ! wget http://ccb.nki.nl/software/aces/ACES.tar.gz
% ! tar -xvzf ACES.tar.gz

%% Initialization
clc;
clear;
h5f_All = './ACES/experiments/data/U133A_combat.h5';
h5f_PerStudy = './ACES/experiments/data/U133A_combat_separate_studies.h5';
xls_Clic = './ACES/experiments/data/@expdesc_breast_3597_qc_110916.xls';

%% Loading full data
fprintf('Reading full dataset from [%s].\n', h5f_All);
Gene_Expression = h5read(h5f_All, '/U133A_combat_RFS/ExpressionData')';
Gene_Entrez = h5read(h5f_All, '/U133A_combat_RFS/GeneLabels');
Patient_Class_Labels = h5read(h5f_All, '/U133A_combat_RFS/PatientClassLabels');
Patient_ID = deblank(h5read(h5f_All, '/U133A_combat_RFS/PatientLabels'));
[n_pat, n_gene] = size(Gene_Expression);
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
	Batch_PatID = deblank(h5read(h5f_PerStudy, [grp_name '/PatientLabels']));
	if ~isequal(Batch_Entz, Gene_Entrez), error(); end
	
	%% Check Expressions
	for pi=1:numel(Batch_PatID)
		pat_ind = Pat_Map(Batch_PatID{pi});
		if ~isequal(Batch_Expr(pi,:), Gene_Expression(pat_ind,:)) || ~isequal(Batch_PatID{pi}, Patient_ID{pat_ind})
			error();
		end
		Study_Index(pat_ind) = step;
	end
	step = step + 1;
end
Gene_Entrez = deblank(Gene_Entrez);

%% Load XSL Clinical
fprintf('Reading clinical data from [%s]\n', xls_Clic);
Clic_tbl = readtable(xls_Clic);
n_item = size(Clic_tbl, 1);
Patient_Info = struct('PatientID', Patient_ID, ...
					  'StudyName', Study_Name(Study_Index), ...
					  'AcesPatientClassLabel', num2cell(strcmp(Patient_Class_Labels, 'TRUE')));
for ti=1:n_item
	if Pat_Map.isKey(Clic_tbl.AffyID{ti})
		pat_ind = Pat_Map(Clic_tbl.AffyID{ti});
		if ~isequal(Patient_Info(pat_ind).PatientID, Clic_tbl.AffyID{ti}), error(); end
		Patient_Info(pat_ind).StudyGSE = Clic_tbl.DataSet{ti};
		Patient_Info(pat_ind).PatientOrgID = Clic_tbl.PatientID{ti};
		Patient_Info(pat_ind).ERStatus = str2double(Clic_tbl.ERStatus{ti});
		Patient_Info(pat_ind).ERStatusArray = Clic_tbl.ER_status_array(ti);
		Patient_Info(pat_ind).ESR1ExpressionOnArray = str2double(Clic_tbl.ESR1ExpressionOnArray(ti));
		Patient_Info(pat_ind).PGRStatus = nan;
		Patient_Info(pat_ind).Her2StatusOnArray = Clic_tbl.HER2_status_array(ti);
		Patient_Info(pat_ind).Her2ExpressionOnArray = Clic_tbl.HER2ExpressionOnArray(ti);
		Patient_Info(pat_ind).TumorSize = Clic_tbl.Size(ti);
		Patient_Info(pat_ind).LymphNodeStatus = str2double(Clic_tbl.LymphNodeStatus{ti});
		Patient_Info(pat_ind).Age = str2double(Clic_tbl.Age{ti});
		Patient_Info(pat_ind).Grade = str2double(Clic_tbl.Grade{ti});
		Patient_Info(pat_ind).DMFSTime = floor(str2double(Clic_tbl.DMFS_time{ti}) * 365);
		Patient_Info(pat_ind).DMFSEvent = str2double(Clic_tbl.DMFS_event_1_relapse_{ti});
		Patient_Info(pat_ind).RFSTime = floor(str2double(Clic_tbl.RFS_time{ti}) * 365);
		Patient_Info(pat_ind).RFSEvent = str2double(Clic_tbl.RFS_event_1_relapse_{ti});
		Patient_Info(pat_ind).OSTime = 0;
		Patient_Info(pat_ind).OSEvent = 0;
		Patient_Info(pat_ind).Tissue = 0;
		Patient_Info(pat_ind).Treatment = Clic_tbl.Treatment{ti};
		Patient_Info(pat_ind).NoTreatment = Clic_tbl.No_treatment(ti);
		Patient_Info(pat_ind).Platform = Clic_tbl.Platform{ti};
		Patient_Info(pat_ind).DeathTime = floor(str2double(Clic_tbl.Death_time{ti}) * 365);
		Patient_Info(pat_ind).DeathEvent = str2double(Clic_tbl.Death_event_1_death_{ti});
		Patient_Info(pat_ind).DeathFromBC = str2double(Clic_tbl.DeathFromBC{ti});
		Patient_Info(pat_ind).SEER_500 = Clic_tbl.SEER_500(ti);
		Patient_Info(pat_ind).EREndocrine = Clic_tbl.ER_endocrine(ti);
		Patient_Info(pat_ind).TripleNegative = Clic_tbl.Triple_negative(ti);
	end
end
% if ~isequal(Patient_ID, {Patient_Info.PatientID}')

%% Run on python
% import cPickle
% from scipy.io import savemat
% fid = open('./ACES/experiments/data/U133A_combat_classlabels.pickle')
% ClinicalFeatures = cPickle.load(fid)
% fid.close()
% savemat('./ACES_ClinicalFeatures.mat', {'ClinicalFeatures':ClinicalFeatures})

%% Load clinical variables
load('./ACES_ClinicalFeatures.mat');
Clic_PatID = deblank(cellstr(ClinicalFeatures{2}));
Clic_Unknown = [ClinicalFeatures{3:6} ClinicalFeatures{7}' ClinicalFeatures{8}{2}];
Clic_Subtype = deblank(cellstr(ClinicalFeatures{8}{1}));
for pi=1:numel(Clic_PatID)
	if Pat_Map.isKey(Clic_PatID{pi})
		pat_ind = Pat_Map(Clic_PatID{pi});
		Patient_Info(pat_ind).CancerSubtype = Clic_Subtype{pi};
		Patient_Info(pat_ind).Unknown = Clic_Unknown(pi,:);
	end
end

%% Saving Data
sav_name = 'ACES_Combined.mat';
fprintf('Saving data in [%s]\n', sav_name);
save(sav_name, 'Gene_Expression', 'Patient_Info', 'Gene_Entrez');

