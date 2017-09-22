%% Download data
% ! wget http://ccb.nki.nl/software/aces/ACES.tar.gz
% ! tar -xvzf ACES.tar.gz

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
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
Clic_Map = containers.Map(Clic_tbl.AffyID, 1:n_item);
Patient_Info = struct('PatientID', Patient_ID, ...
					  'StudyName', Study_Name(Study_Index), ...
					  'AcesPatientClassLabel', num2cell(strcmp(Patient_Class_Labels, 'TRUE')));
for pi=1:n_pat
	if Clic_Map.isKey(Patient_ID{pi})
		pat_ind = Clic_Map(Patient_ID{pi});
		if ~isequal(Patient_Info(pi).PatientID, Clic_tbl.AffyID{pat_ind}), error(); end
		Patient_Info(pi).StudyGSE = Clic_tbl.DataSet{pat_ind};
		Patient_Info(pi).PatientOrgID = Clic_tbl.PatientID{pat_ind};
		Patient_Info(pi).ERStatus = str2double(Clic_tbl.ERStatus{pat_ind});
		Patient_Info(pi).ERStatusArray = Clic_tbl.ER_status_array(pat_ind);
		Patient_Info(pi).ESR1ExpressionOnArray = str2double(Clic_tbl.ESR1ExpressionOnArray(pat_ind));
		Patient_Info(pi).PGRStatus = nan;
		Patient_Info(pi).Her2StatusOnArray = Clic_tbl.HER2_status_array(pat_ind);
		Patient_Info(pi).Her2ExpressionOnArray = Clic_tbl.HER2ExpressionOnArray(pat_ind);
		Patient_Info(pi).TumorSize = Clic_tbl.Size(pat_ind);
		Patient_Info(pi).LymphNodeStatus = str2double(Clic_tbl.LymphNodeStatus{pat_ind});
		Patient_Info(pi).Age = str2double(Clic_tbl.Age{pat_ind});
		Patient_Info(pi).Grade = str2double(Clic_tbl.Grade{pat_ind});
		Patient_Info(pi).DMFSTime = floor(str2double(Clic_tbl.DMFS_time{pat_ind}) * 365);
		Patient_Info(pi).DMFSEvent = str2double(Clic_tbl.DMFS_event_1_relapse_{pat_ind});
		Patient_Info(pi).RFSTime = floor(str2double(Clic_tbl.RFS_time{pat_ind}) * 365);
		Patient_Info(pi).RFSEvent = str2double(Clic_tbl.RFS_event_1_relapse_{pat_ind});
		Patient_Info(pi).OSTime = nan;
		Patient_Info(pi).OSEvent = nan;
		Patient_Info(pi).Tissue = '';
		Patient_Info(pi).Treatment = Clic_tbl.Treatment{pat_ind};
		Patient_Info(pi).NoTreatment = Clic_tbl.No_treatment(pat_ind);
		Patient_Info(pi).Platform = Clic_tbl.Platform{pat_ind};
		Patient_Info(pi).DeathTime = floor(str2double(Clic_tbl.Death_time{pat_ind}) * 365);
		Patient_Info(pi).DeathEvent = str2double(Clic_tbl.Death_event_1_death_{pat_ind});
		Patient_Info(pi).DeathFromBC = str2double(Clic_tbl.DeathFromBC{pat_ind});
		Patient_Info(pi).SEER_500 = Clic_tbl.SEER_500(pat_ind);
		Patient_Info(pi).EREndocrine = Clic_tbl.ER_endocrine(pat_ind);
		Patient_Info(pi).TripleNegative = Clic_tbl.Triple_negative(pat_ind);
	else
		%% Filling up the empty items
		for fn=fieldnames(Patient_Info(pi))'
			if numel(Patient_Info(pi).(fn{1}))==0
				if isnumeric(Patient_Info(1).(fn{1}))
					Patient_Info(pi).(fn{1}) = nan;
				else
					Patient_Info(pi).(fn{1}) = '';
				end
			end
		end
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
Patient_Info = struct2table(Patient_Info);

%% Saving Data
sav_name = 'ACES_Combined.mat';
fprintf('Saving data in [%s]\n', sav_name);
save(sav_name, 'Gene_Expression', 'Patient_Info', 'Gene_Entrez');

