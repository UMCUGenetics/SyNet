
%% Initialization
clc;
clear;
data_lst = {'ACES' 'HaibeKains' 'METABRIC' 'TCGA'};
Header_List = {
%%	SyNet					ACES							HAIBE				METABRIC				TCGA
	'SurvivalTime'			''								''					''						''
	'Prognostic_Status'		'AcesPatientClassLabel'			''					''						''
	'Age'					'Age'							'Age'				'age_at_diagnosis'		'Age_at_Initial_Pathologic_Diagnosis_nature2012'
	'Subtype'				'CancerSubtype'					'Subtype'			'Pam50Subtype'			'Subtype'
	'ERStatus'				'ERStatus'						'ERStatus'			'ER_Expr'				'ER_Status_nature2012'
	'Her2Status'			'Her2StatusOnArray'				'Her2Status'		'Her2_Expr'				'HER2_Final_Status_nature2012'
	'PGRStatus'				'PGRStatus'						'PGRStatus'			'PR_Expr'				'PR_Status_nature2012'
	'Grade'					'Grade'							'Grade'				'grade'					''
	'LymphNodeStatus'		'LymphNodeStatus'				'LymphNodeStatus'	'lymph_nodes_positive'	'Node_Coded_nature2012'
	'TumorSize'				'TumorSize'						'TumorSize'			'size'					''	
	'Platform'				'Platform'						'Platform'			'Platform'				'Platform'
	'RFSEvent'				'RFSEvent'						'RFSEvent'			''						''
	'RFSTime'				'RFSTime'						'RFSTime'			''						'x_RFS'
	'DMFSEvent'				'DMFSEvent'						'DMFSEvent'			'dss_event'				''
	'DMFSTime'				'DMFSTime'						'DMFSTime'			''						''
	'OSEvent'				'OSEvent'						'OSEvent'			'os_event'				'OS_event_nature2012'
	'OSTime'				'OSTime'						'OSTime'			'survival_time'			'OS_Time_nature2012'
	'PatientID'				'PatientID'						'PatientID'			'patient_id'			'sampleID'
	'Treatment'				'Treatment'						'Treatment'			'Treatment'				'history_of_neoadjuvant_treatment'
	'StudyName'				'StudyName'						'StudyName'			''						''
	'StudyGSE'				'StudyGSE'						''					''						''
};

%% Loading ACES data
fprintf('Reading ACES data.\n');
data_aces = load('../ACES/ACES_Combined.mat');
data_aces.Patient_Info = SelectFromTable(data_aces.Patient_Info, Header_List(:,2), Header_List(:,1));
data_aces.Patient_Info.SurvivalTime = data_aces.Patient_Info.RFSTime;
data_aces.Gene_Entrez = strrep(data_aces.Gene_Entrez, 'Entrez_', '');

%% Loading HaibeKains data
fprintf('Reading Haibe data.\n');
data_haib = load('../HaibeKains/HaibeKains_Combined.mat');
is_dup = data_haib.Patient_Info.IsDuplicated == 1;
data_haib.Patient_Info(is_dup, :) = [];
data_haib.Gene_Expression(is_dup, :) = [];
data_haib.Patient_Info = SelectFromTable(data_haib.Patient_Info, Header_List(:,3), Header_List(:,1));
data_haib.Patient_Info = CastFields2Num(data_haib.Patient_Info, Header_List([3 5:10 12:17],1));
data_haib.Patient_Info = getSurvivalTime(data_haib.Patient_Info);
data_haib.Patient_Info.Prognostic_Status = data_haib.Patient_Info.SurvivalTime <= 1825;

%% Loading METABRIC data
fprintf('Reading Metabric data.\n');
data_meta = load('../METABRIC/METABRIC_Combined.mat');
data_meta.Patient_Info = SelectFromTable(data_meta.Patient_Info, Header_List(:,4), Header_List(:,1));
data_meta.Patient_Info = RepFieldsWithValue(data_meta.Patient_Info, Header_List(5:7,1), {'\+' '\-'}, {'1' '0'});
data_meta.Patient_Info = CastFields2Num(data_meta.Patient_Info, Header_List(5:7,1));
data_meta.Patient_Info.SurvivalTime = data_meta.Patient_Info.OSTime;
data_meta.Patient_Info.Prognostic_Status = data_meta.Patient_Info.SurvivalTime <= 1825;
data_meta.Patient_Info.StudyName = repmat({'METABRIC'}, size(data_meta.Patient_Info,1), 1);
data_meta.Patient_Info.StudyGSE = repmat({'NA'}, size(data_meta.Patient_Info,1), 1);

%% Converting Entrez to Hugo GeneName
% ! wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
% # For complete description refer to ftp://ftp.ncbi.nih.gov/gene/DATA/README
% ! zcat Homo_sapiens.gene_info | awk '{print($2"\t"$3"\t"$5)}' > Entrez2GeneName.txt
fid = fopen('./Entrez2GeneName.txt', 'r');
f_cell = textscan(fid, '%s%s%s', 'Delimiter', '\t', 'HeaderLines', 1, 'ReturnOnError', 0);
if ~feof(fid), error(); end
Ent2Name = containers.Map(f_cell{1}, f_cell{2});
Syn2Ent = containers.Map(f_cell{2}, f_cell{1});
Gene_Synm = regexp(f_cell{3}, '\|', 'split');
for si=1:numel(f_cell{3})
	gene_lst = Gene_Synm{si};
	for gi=1:numel(gene_lst)
		if ~strcmp(gene_lst{gi}, '-') && ~Syn2Ent.isKey(gene_lst{gi}) % If gene has an Entrez already, we keep it
			Syn2Ent(gene_lst{gi}) = f_cell{1}{si};
		end
	end
end

%% Load TCGA data
fprintf('Reading TCGA data.\n');
data_tcga = load('../TCGA/TCGA_Combined.mat');
data_tcga.Patient_Info = SelectFromTable(data_tcga.Patient_Info, Header_List(:,5), Header_List(:,1));
data_tcga.Patient_Info = RepFieldsWithValue(data_tcga.Patient_Info, Header_List([5:7 9],1), {'Positive' 'Negative' 'Equivocal'}, {'1' '0' 'NA'});
data_tcga.Patient_Info = CastFields2Num(data_tcga.Patient_Info, Header_List([5:7 9],1));
% data_tcga.Patient_Info.Age = str2double(data_tcga.Patient_Info.Age);
data_tcga.Patient_Info = getSurvivalTime(data_tcga.Patient_Info);
is_normal = ~cellfun('isempty', regexp(data_tcga.Patient_Info.PatientID, '-11$'));
data_tcga.Patient_Info(is_normal, :) = [];
data_tcga.Gene_Expression(is_normal, :) = [];
data_tcga.Patient_Info.StudyName = repmat({'TCGA'}, size(data_tcga.Patient_Info,1), 1);
data_tcga.Patient_Info.StudyGSE = repmat({'NA'}, size(data_tcga.Patient_Info,1), 1);
n_gene = numel(data_tcga.Gene_Name);
data_tcga.Gene_Entrez = cell(n_gene, 1);
for gi=1:n_gene
	if Syn2Ent.isKey(data_tcga.Gene_Name{gi})
		data_tcga.Gene_Entrez{gi} = Syn2Ent(data_tcga.Gene_Name{gi});
	else
		data_tcga.Gene_Entrez{gi} = '';
	end
end

%% Unify expression
data_syne.Gene_Entrez = intersect(data_aces.Gene_Entrez, data_haib.Gene_Entrez);
data_syne.Gene_Entrez = intersect(data_syne.Gene_Entrez, data_meta.Gene_Entrez);
data_syne.Gene_Entrez = intersect(data_syne.Gene_Entrez, data_tcga.Gene_Entrez);
Expr_aces = UnifyExpression(data_aces, data_syne.Gene_Entrez);
Expr_haib = UnifyExpression(data_haib, data_syne.Gene_Entrez);
Expr_meta = UnifyExpression(data_meta, data_syne.Gene_Entrez);
Expr_tcga = UnifyExpression(data_tcga, data_syne.Gene_Entrez);
data_syne.Gene_Expression = [Expr_aces; Expr_haib; Expr_meta; Expr_tcga];
data_syne.Patient_Info = [data_aces.Patient_Info; data_haib.Patient_Info; data_meta.Patient_Info; data_tcga.Patient_Info];

%% Loading clinical data
Item_Info = readtable(fn_clic, 'HeaderLines', 0, 'TreatAsEmpty', 'NA');
n_item = size(Item_Info.sampleID, 1);
Item_Map = containers.Map(Item_Info.sampleID, 1:n_item);
Patient_Info = table();
for pi=1:n_pat
	item_ind = Item_Map(Patient_ID{pi});
	Patient_Info(pi, :) = Item_Info(item_ind, :);
end
if ~isequal(Patient_ID, Patient_Info.sampleID)
	error();
end

%% Replacing nans with median
for gi=1:n_gene
	is_nan = isnan(Gene_Expression(:, gi));
	if any(is_nan)
		Gene_Expression(is_nan, gi) = median(Gene_Expression(~is_nan, gi));
	end
end
if any(isnan(Gene_Expression(:))), error(); end

%% Saving Data
sav_name = 'TCGA_Combined.mat';
fprintf('Saving data in [%s]\n', sav_name);
save(sav_name, 'Gene_Expression', 'Patient_Info', 'Gene_Name');

%% ///////////////////// Functions
function out_tbl = SelectFromTable(in_tbl, SrcHeader, TarHeader)
n_row = size(in_tbl, 1);
for hi=1:numel(SrcHeader)
	if strcmp(SrcHeader{hi}, '')
		out_tbl(:, TarHeader{hi}) = table(nan(n_row, 1));
	else
		out_tbl(:, TarHeader{hi}) = in_tbl(:, SrcHeader{hi});
	end
end
end

function tbl = getSurvivalTime(tbl)
for ti=1:size(tbl, 1)
	if ~isnan(tbl.RFSTime(ti))
		tbl.SurvivalTime(ti) = tbl.RFSTime(ti);
	elseif ~isnan(tbl.DMFSTime(ti))
		tbl.SurvivalTime(ti) = tbl.DMFSTime(ti);
	else
		tbl.SurvivalTime(ti) = tbl.OSTime(ti);
	end
end
end

function Gene_Expression = UnifyExpression(tbl, Gene_Entrez)
n_targene = numel(Gene_Entrez);
[n_pat, n_srcgene] = size(tbl.Gene_Expression);
Gene_Map = containers.Map(tbl.Gene_Entrez, 1:n_srcgene);
Gene_Expression = nan(n_pat, n_targene);
for gi=1:n_targene
	if Gene_Map.isKey(Gene_Entrez{gi})
		gene_ind = Gene_Map(Gene_Entrez{gi});
		Gene_Expression(:, gi) = tbl.Gene_Expression(:, gene_ind);
	end
end
end

function tbl = CastFields2Num(tbl, field_lst)
fprintf('Casting fields to number: ');
for fi=1:numel(field_lst)
	fprintf('%s, ', field_lst{fi});
	tbl.(field_lst{fi}) = str2double(tbl.(field_lst{fi}));
end
fprintf('\n');
end

function tbl = RepFieldsWithValue(tbl, field_lst, src_str, tar_str)
fprintf('Replacing fields values: ');
for fi=1:numel(field_lst)
	fprintf('%s, ', field_lst{fi});
	for si=1:numel(src_str)
		tbl.(field_lst{fi}) = regexprep(tbl.(field_lst{fi}), ['^' src_str{si} '$'], tar_str{si});
	end
end
end

function tbl = ReplaceNAN(tbl, field_lst, str)
fprintf('Replacing fields [%s]: ', str{1});
for fi=1:numel(field_lst)
	fprintf('%s, ', field_lst{fi});
	if any(~isnan(tbl.(field_lst{fi})))
		error();
	end
	tbl.(field_lst{fi}) = repmat(str, size(tbl,1), 1);
end
end
