clc;
clear

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
data_path = './csv_data/';

%% Load study info
Study_Info = readtable('./csv_data/Data_DDB_ddb.csv', 'HeaderLines', 0);
Study_Info.Properties.VariableNames = {'RowName' 'StudyName' 'Citation' 'Platform' 'Reference' 'DFSTime' 'OSTime' 'Auto'};
if ismac, Study_Info = Study_Info(1:3, :); end
n_study = size(Study_Info, 1);

%% Load data files
Batch_Pat = cell(n_study, 1);
Batch_Prb = cell(n_study, 1);
Batch_Expr = cell(n_study, 1);
Gene_Entrez = cell(n_study, 1);
for si=1:n_study
	fprintf('[%d/%d] Reading [%s] study:\n', si, n_study, Study_Info.StudyName{si});
	
	%% Read files
	Batch_Pat{si} = readClinical([data_path 'Data_' Study_Info.StudyName{si} '_demo.csv'], Study_Info(si,:));
	Batch_Prb{si} = readProb([data_path 'Data_' Study_Info.StudyName{si} '_annot.csv']);
	[Batch_Expr{si}, Batch_PatID, BatchPrbID] = readExpression([data_path 'Data_' Study_Info.StudyName{si} '_data.csv']);
	if ~isequal({Batch_Pat{si}.RowName}', Batch_PatID) || ~isequal({Batch_Prb{si}.RowName}', BatchPrbID)
		error('Data are not consistent');
	end
	
	%% Save Probs
	Gene_Entrez{si} = {Batch_Prb{si}.EntrezID};
	fprintf('=====\n');
end
fprintf('\n');

%% Map probs to genes
fprintf('Mapping probes to genes ...\n');
Gene_Entrez = unique([Gene_Entrez{:}])';
Gene_Entrez(strcmp(Gene_Entrez, 'NA')) = [];
n_Entz = numel(Gene_Entrez);
GMap = containers.Map(Gene_Entrez, 1:n_Entz);
Gene_Expression = zeros(0, n_Entz);
Patient_Info = [];
Prob_ID = cell(n_Entz, 1);
Gene_Name = cell(n_Entz, 1);
for si=1:n_study
	fprintf('[%d/%d] Adjusting [%s] study:\n', si, n_study, Study_Info.StudyName{si});
	Batch_Entrz = {Batch_Prb{si}.EntrezID}';
	Prb_grp = accumarray(grp2idx(Batch_Entrz), 1:numel(Batch_Entrz), [], @(i) {i});
	n_grp = numel(Prb_grp);
	n_pat = size(Batch_Expr{si}, 1);
	tmp_Expr = nan(n_pat, n_Entz);
	for gi=1:n_grp
		showprogress(gi, n_grp);
		if GMap.isKey(Batch_Prb{si}(Prb_grp{gi}(1)).EntrezID)
			%% Check for nan ratio
			prb_ind = Prb_grp{gi};
			for pi=1:numel(prb_ind)
				is_nan = isnan(Batch_Expr{si}(:, prb_ind(pi)));
				if sum(is_nan)/n_pat>0.75
					prb_ind(pi) = nan;
				else
					Batch_Expr{si}(is_nan, prb_ind(pi)) = median(Batch_Expr{si}(~is_nan, prb_ind(pi)));
				end
			end
			if sum(~isnan(prb_ind))<1
				fprintf('Warning: No probs left for [%s] Entrez, [%s] gene. [%0.1f%%] are NANs\n', ...
					Batch_Prb{si}(Prb_grp{gi}(1)).EntrezID, Batch_Prb{si}(Prb_grp{gi}(1)).GeneName, sum(is_nan)*100/n_pat);
				continue;
			end
			Prb_grp{gi} = prb_ind(~isnan(prb_ind));
			
			%% Get the top variable prob
			prob_set = Batch_Prb{si}(Prb_grp{gi});
			prb_std = std(Batch_Expr{si}(:, Prb_grp{gi}));
			if any(isnan(prb_std)), error('Std of a prob is NaN.'); end
			[~, sid] = sort(prb_std, 'Descend');
			prob_set = prob_set(sid);
			Prb_grp{gi} = Prb_grp{gi}(sid);
			
			ge_ind = GMap(prob_set(1).EntrezID);
			pr_ind = Prb_grp{gi}(1);
			tmp_Expr(:, ge_ind) = Batch_Expr{si}(:, pr_ind);
			Prob_ID{ge_ind} = [Prob_ID{ge_ind} {prob_set.ProbID}];
			Gene_Name{ge_ind} = [Gene_Name{ge_ind} {prob_set.GeneName}];
		end
	end
	
	%% Add to collection
	Gene_Expression = [Gene_Expression; tmp_Expr];
	Patient_Info = [Patient_Info; Batch_Pat{si}];
end

%% Combine probs
fprintf('Combining prob IDs ...\n');
for gi=1:n_Entz
	if numel(Prob_ID{gi})<1
		fprintf('Warning: No probs is used for [%s] entrez.\n', Gene_Entrez{gi});
	else
		Prob_ID{gi} = strjoin(unique(Prob_ID{gi}), ';');
		Gene_Name{gi} = strjoin(unique(Gene_Name{gi}), ';');
	end
end

%% Saving Data
sav_name = 'HaibeKains_Combined.mat';
fprintf('Saving data in [%s]\n', sav_name);
save(sav_name, 'Gene_Expression', 'Patient_Info', 'Prob_ID', 'Gene_Name', 'Gene_Entrez', 'Study_Info');

%% Functions \\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
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
ign_lst = {};
for hi=1:n_Header
	showprogress(hi, n_Header);
	ind = find(ismember(Struct_Headers(:,1), f_Header{hi}));
	if numel(ind)==0
		ign_lst = [ign_lst f_Header{hi}];
		continue;
	elseif numel(ind)>1
		error('Duplicated header found.\n');
	end
	
	[Prob_Info.(Struct_Headers{ind,2})] = deal(f_cell{hi}{:});
end
if numel(ign_lst)>0
	fprintf('Warning: [%s] fields are ignored.\n', strjoin(ign_lst, ', '));
end
fprintf('[%0.1f%%] of probs have EntrezID.\n', numel(unique({Prob_Info(:).EntrezID}))*100/numel(Prob_Info));
end

function Patient_Info = readClinical(csv_name, Study_Info)
Struct_Headers = {
	'samplename'	'PatientID' %1
	'dataset'		'StudyName' %2
	'series'		'Series' %3
	'id'			'PatientOrgID' %4
	'er'			'ERStatus' %5
	'pgr'			'PGRStatus' %6
	'her2'			'Her2Status' %7
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
[Patient_Info(:).Platform] = deal(Study_Info.Platform);
end
