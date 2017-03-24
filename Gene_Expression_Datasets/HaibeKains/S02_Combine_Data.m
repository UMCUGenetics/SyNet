clc;
clear

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
data_path = './csv_data/';
data_lst = {'UPP' 'VDX'};
n_data = numel(data_lst);

%% Load data files
Batch_Pat = cell(n_data, 1);
Batch_Prb = cell(n_data, 1);
Batch_Expr = cell(n_data, 1);
All_Entrez = cell(n_data, 1);
for di=1:n_data
	fprintf('[%d/%d] Reading [%s] study:\n', di, n_data, data_lst{di});
	
	%% Read files
	Batch_Pat{di} = readClinical([data_path 'Data_' data_lst{di} '_demo.csv']);
	Batch_Prb{di} = readProb([data_path 'Data_' data_lst{di} '_annot.csv']);
	[Batch_Expr{di}, Batch_PatID, BatchPrbID] = readExpression([data_path 'Data_' data_lst{di} '_data.csv']);
	if ~isequal({Batch_Pat{di}.RowName}', Batch_PatID) || ~isequal({Batch_Prb{di}.RowName}', BatchPrbID)
		error('Data are not consistent');
	end
	
	%% Save Probs
	All_Entrez{di} = {Batch_Prb{di}.EntrezID};
end

%% Map probs to genes
fprintf('Mapping probes to genes ...\n');
All_Entrez = unique([All_Entrez{:}])';
All_Entrez(strcmp(All_Entrez, 'NA')) = [];
n_All_Genes = numel(All_Entrez);
GMap = containers.Map(All_Entrez, 1:n_All_Genes);
GeneExpression = zeros(0, n_All_Genes);
Patient_Info = [];
Prob_ID = cell(n_All_Genes, 1);
for di=1:n_data
	fprintf('[%d/%d] Adjusting [%s] study:\n', di, n_data, data_lst{di});
	Batch_Entrz = {Batch_Prb{di}.EntrezID}';
	[n_pat, n_gene] = size(Batch_Expr{di});
	Prb_grp = accumarray(grp2idx(Batch_Entrz), 1:n_gene, [], @(i) {i});
	n_grp = numel(Prb_grp);
	tmp_Expr = nan(n_pat, n_All_Genes);
	for gi=1:n_grp
		showprogress(gi, n_grp);
		prob_set = Batch_Prb{di}(Prb_grp{gi});
		if GMap.isKey(prob_set(1).EntrezID)
			[~, sid] = sort(std(Batch_Expr{di}(:, Prb_grp{gi})), 'Descend');
			prob_set = prob_set(sid);
			Prb_grp{gi} = Prb_grp{gi}(sid);
			
			ge_ind = GMap(prob_set(1).EntrezID);
			pr_ind = Prb_grp{gi}(1);
			tmp_Expr(:, ge_ind) = Batch_Expr{di}(:, pr_ind);
		end
	end
end

%% Order Genes
fprintf('Ordering the data.');
for di=1:n_data
	showprogress(di, n_data);
	
	[n_pat, n_gene] = size(Batch_Expr{di});
	Gene_Ind = zeros(n_gene, 1);
	for gi=1:n_gene
		Gene_Ind(gi) = GMap(Batch_Prb{di}(gi).GeneName);
	end
	Bat_Exp = zeros(n_pat, GMap.Count);
	Bat_Exp(:) = Batch_Expr{di}(:, Gene_Ind);
end
if isempty(PatientID)
	PatientInfo = Batch_Pat;
	ProbInfo = Batch_Prb;
	ExpressionData = Batch_Expr;
	
	for pi=1:numel(Batch_PatID)
		s
		if isempty(ProbInfo)
			ProbInfo = PI;
		else
			for pi=1:size(PI,1)
			end
		end
	end
end

function [ExpressionData, PatientID, ProbID] = readExpression(csv_name)
fprintf('Reading expressions from [%s]: ', csv_name);
fid = fopen(csv_name, 'r');
ProbID = regexp(fgetl(fid), '\t', 'split')';
n_Prob = numel(ProbID);
f_cell = textscan(fid, ['%s' repmat('%f', 1, n_Prob)], 'HeaderLines', 0, 'Delimiter', '\t', 'TreatAsEmpty', {'NA'}, 'ReturnOnError', 0);
if ~feof(fid), error('File is corrupted ...'); end
fclose(fid);
n_Patient = size(f_cell{1}, 1);

PatientID = f_cell{1};
f_cell(1) = [];
ExpressionData = zeros(n_Patient, n_Prob);
for pi=1:n_Prob
	showprogress(pi, n_Prob);
	ExpressionData(:, pi) = f_cell{pi};
end
end

function Prob_Info = readProb(csv_name)
Struct_Headers = {
	'probe'					'ProbID' %1
	'Gene.title'			'GeneFullName' %2
	'Gene.symbol'			'GeneName' %3
	'Gene.ID'				'EntrezIDList' %4
	'EntrezGene.ID'			'EntrezID' %5
	'UniGene.title'			'UniGeneTitle' %6
	'UniGene.symbol'		'UniGeneSymbol' %7
	'UniGene.ID'			'UniGeneID' %8
	'Nucleotide.Title'		'NucleotideTitle' %9
	'GI'					'GI' %10
	'GenBank.Accession'		'GenBankAccession' %11
	'Platform_CLONEID'		'PlatformCLONEID' %12
	'Platform_ORF'			'PlatformORF' %13
	'Platform_SPOTID'		'Platform_SPOTID' %14
	'Chromosome.location'	'ChromosomeLocation' %15
	'Chromosome.annotation'	'ChromosomeAnnotation' %16
	'GO.Function'			'GOFunction' %17
	'GO.Process'			'GOProcess' %18
	'GO.Component'			'GOComponent' %19
	'GO.Function.1'			'GOFunction1' %20
	'GO.Process.1'			'GOProcess1' %21
	'GO.Component.1'		'GOComponent1' %22
	};

fprintf('Reading prob information from [%s]: ', csv_name);
fid = fopen(csv_name, 'r');
f_Header = regexp(fgetl(fid), '\t', 'split')';
n_Header = numel(f_Header);
f_cell = textscan(fid, repmat('%s', 1, n_Header+1), 'HeaderLines', 0, 'Delimiter', '\t', 'TreatAsEmpty', {'NA'}, 'ReturnOnError', 0);
if ~feof(fid), error('File is corrupted ...'); end
fclose(fid);
n_Prob = size(f_cell{1}, 1);

Prob_Info = repmat(struct(), n_Prob, 1);
[Prob_Info(:).RowName] = deal(f_cell{1}{:});
f_cell(1) = [];
for hi=1:n_Header
	showprogress(hi, n_Header);
	ind = find(ismember(Struct_Headers(:,1), f_Header{hi}));
	if numel(ind)~=1, error(); end
	
	[Prob_Info.(Struct_Headers{ind,2})] = deal(f_cell{hi}{:});
end
end

function Patient_Info = readClinical(csv_name)
Struct_Headers = {
	'samplename'	'PatientID' %1
	'dataset'		'StudyName' %2
	'series'		'Series' %3
	'id'			'PatientOrgID' %4
	'er'			'ERStatus' %5
	'pgr'			'PGRStatus' %6
	'her2'			'HerTStatus' %7
	'size'			'TumorSize' %8
	'node'			'LymphNodeStatus' %9
	'age'			'Age' %10
	'grade'			'Grade' %11
	't.dmfs'		'DMFSTime' %12
	'e.dmfs'		'DMFSEvent' %13
	't.rfs'			'RFSTime' %14
	'e.rfs'			'RFSEvent' %15
	't.os'			'OSTime' %16
	'e.os'			'OSEvent' %17
	'tissue'		'Tissue' %18
	'treatment'		'Treatment' %19
	};

fprintf('Reading clinical variables from [%s]: ', csv_name);
fid = fopen(csv_name, 'r');
f_Header = regexp(fgetl(fid), '\t', 'split')';
n_Header = numel(f_Header);
f_cell = textscan(fid, repmat('%s', 1, n_Header+1), 'HeaderLines', 0, 'Delimiter', '\t', 'TreatAsEmpty', {'NA'}, 'ReturnOnError', 0);
n_Patient = numel(f_cell{1});

Patient_Info = repmat(struct(), n_Patient, 1);
[Patient_Info(:).RowName] = deal(f_cell{1}{:});
f_cell(1) = [];
for hi=1:n_Header
	showprogress(hi, n_Header);
	ind = find(ismember(Struct_Headers(:,1), f_Header{hi}));
	if numel(ind)~=1, error(); end
	
	[Patient_Info.(Struct_Headers{ind,2})] = deal(f_cell{hi}{:});
end
end
