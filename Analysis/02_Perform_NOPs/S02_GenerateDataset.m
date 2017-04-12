function ds_id = S02_GenerateDataset(cv_id, net_name)
%% ####
% cv_id = '1608231618';
% net_name = 'STRING';
% net_name = 'DSN-ACES-99.9';
if ismac
	fprintf('*** Warning!: Running on debug mode.\n');
	cv_id = '170411001004';
	net_name = 'DSN-SyNetS1-T20';
end

%% Initialization
%clc;
addpath('../../../../Useful_Sample_Codes/ShowProgress');
cv_path = './CV_Files/';
dataset_path = './Dataset_Files/';
%rng shuffle % creates a different seed each time, but its not nessasary as seed is set in the job submission
if ~exist('ds_id','var'); ds_id=[datestr(now, 'yymmddhhMM') num2str(randi(999999), '%06d')]; end

%% Get CV info
cv_list = dir([cv_path 'CID-' cv_id '_*.mat']);
if numel(cv_list)~=1, error('Missing or duplicated CV info found. [%s]', strjoin({cv_list.name}, ', ')); end
cv_name = [cv_path cv_list.name];
fprintf('Loading CV info [%s] ...\n', cv_name);
load(cv_name);
% if ~isempty(strfind(net_name, te_info.GEName)), error('Test data can not be used in training network.\n'); end

dataset_name = sprintf([dataset_path cv_list.name(1:end-4) '_NET-' net_name '_DID-%s.mat'], ds_id);
if exist(dataset_name, 'file'), error('Dataset is already generated. [%s]', dataset_name); end
fprintf('Generating a dataset in [%s].\n', dataset_name);

%% Get Gene Expression
data = load(tr_info.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
if ~isequal(data.Gene_Name, tr_info.Gene_Name), error(); end
tr_info.Gene_Expression = data.Gene_Expression;
tr_info.Patient_Label = (data.Patient_Label==1)*2-1;

data = load(te_info.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
if ~isequal(data.Gene_Name, te_info.Gene_Name), error(); end
te_info.Gene_Expression = data.Gene_Expression;
te_info.Patient_Label = (data.Patient_Label==1)*2-1;
clear data

%% Get Network Info
net_info = getNetworkInfo(net_name, tr_info);
Valid_Gene_List = intersect(tr_info.Gene_Name, te_info.Gene_Name);
Valid_Gene_List = intersect(Valid_Gene_List, net_info.Gene_Name);

%% Generate dataset
fprintf('Generating training set from [%s, %s] ...\n', tr_info.GEName, net_info.net_name);
DatasetTr = getDataset(tr_info, net_info, Valid_Gene_List);
DatasetTr.iCvPar = tr_info.iCvPar;
DatasetTr.GEName = tr_info.GEName;
DatasetTr.NetName = net_name;

fprintf('Generating test set from [%s] ...\n', te_info.GEName);
DatasetTe = getDataset(te_info, [], Valid_Gene_List);
DatasetTe.GEName = te_info.GEName;

%% Save dataset
fprintf('Saving dataset in [%s].\n', dataset_name);
save(dataset_name, 'DatasetTr', 'DatasetTe');
end

function net_info = getNetworkInfo(net_name, tr_data)
net_info.net_name = net_name;
switch net_name
	case 'Random'
		net_info.Gene_Name = tr_data.Gene_Name;
	case 'Corr'
		net_info.Gene_Name = tr_data.Gene_Name;
	case {'BioGRID' 'CPDB' 'FunCoup-Complex' 'FunCoup-DOM' 'FunCoup-GIN' 'FunCoup-Metabolic' 'FunCoup-MEX' 'FunCoup-MIR' 'FunCoup-PEX' 'FunCoup-PFC' 'FunCoup-PHP' 'FunCoup-PPI' 'FunCoup-SCL' 'FunCoup-Signaling' 'FunCoup-TFB' 'HiC-ES' 'HiC-Full' 'HIN' 'HINT' 'HPRD' 'I2D' 'INstruct' 'iRefIndex' 'KEGG' 'MSigDB' 'Multinet' 'Reactome' 'String-CoEx' 'String-CoOc' 'String-Data' 'String-Exp' 'String-Fus' 'String-Neig' 'String-Txt' 'String' 'UniHI'}
		switch net_name
			case 'BioGRID'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/BioGRID_Full_GN.csv';
			case 'CPDB'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/CPDB_GN_1percent.csv';
			case 'FunCoup-Complex'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_Complex_GN.txt';
			case 'FunCoup-DOM'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_DOM_GN.txt';
			case 'FunCoup-GIN'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_GIN_GN.txt';
			case 'FunCoup-Metabolic'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_Metabolic_GN.txt';
			case 'FunCoup-MEX'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_MEX_GN.txt';
			case 'FunCoup-MIR'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_MIR_GN.txt';
			case 'FunCoup-PEX'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_PEX_GN.txt';
			case 'FunCoup-PFC'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_PFC_GN.txt';
			case 'FunCoup-PHP'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_PHP_GN.txt';
			case 'FunCoup-PPI'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_PPI_GN.txt';
			case 'FunCoup-SCL'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_SCL_GN.txt';
			case 'FunCoup-Signaling'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_Signaling_GN.txt';
			case 'FunCoup-TFB'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/FunCoup_TFB_GN.txt';
			case 'HiC-ES'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/HiC_Full_ES-all-hg19-corrected-200k.tsv';
			case 'HiC-Full'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/HiC_Full.tsv';
			case 'HIN'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/HIN_Full_GN.csv';
			case 'HINT'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/HINT_GN.tsv';
			case 'HPRD'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/HPRD_GN.txt';
			case 'I2D'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/I2D_GN.csv';
			case 'INstruct'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/INstruct_GN.txt';
			case 'iRefIndex'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/iRefIndex_GN.tsv';
			case 'KEGG'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/KEGG_GN.tsv';
			case 'MSigDB'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/MSigDB_GN_1percent.csv';
			case 'Multinet'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/Multinet_GN.txt';
			case 'Reactome'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/Reactome_GN.tsv';
			case 'String-CoEx'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/STRING_coexpression_GN.txt';
			case 'String-CoOc'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/STRING_cooccurence_GN.txt';
			case 'String-Data'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/STRING_database_GN.txt';
			case 'String-Exp'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/STRING_experimental_GN.txt';
			case 'String-Fus'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/STRING_fusion_GN.txt';
			case 'String-Neig'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/STRING_neighborhood_GN.txt';
			case 'String-Txt'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/STRING_textmining_GN.txt';
			case 'String'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/String_Homo_GN_le900t_WeightsRemoved.csv';
			case 'UniHI'
				net_path = '../../../105_Scale_Aware_Topological_Measure_Calculator_GRID_Ready/Graphs/UniHI_GN.txt';
		end
		fid = fopen(net_path, 'r');
		net_cell = textscan(fid, '%s%s', 'Delimiter', '\t', 'ReturnOnError', 0);
		if ~feof(fid), error(); end
		fclose(fid);
		net_info.Gene_Name = unique(vertcat(net_cell{:}));
		net_info.net_cell = net_cell;
		net_info.net_path = net_path;
	otherwise
		if ~isequal(net_name(1:4), 'DSN-') && ~isequal(net_name(1:4), 'CPN-')
			error('Unknown network.');
		end
		net_part = regexp(net_name, '-', 'split');
		net_info.net_name = net_part{1};
		net_info.net_src = net_part{2};
		net_info.param_type = net_part{3}(1);
		net_info.param_val = str2double(net_part{3}(2:end));
		fprintf('Detected [%s] network from [%s], Type: [%s], Value: [%g]\n', net_info.net_name, net_info.net_src, net_info.param_type, net_info.param_val);
		net_path = ['../01_Pairwise_Evaluation_of_Genes/Network_Files/' net_info.net_name '_' net_info.net_src '.mat'];
		load(net_path, 'Net_Adj', 'Gene_Name');
		if ~issymmetric(Net_Adj), error('Adj Matrix is not symetric.\n'); end
		switch net_info.param_type
			case 'T' % Selecting top %d interactions
				scr_val = sort(Net_Adj(:), 'descend');
				adj_tresh = scr_val(net_info.param_val);
			case 'P' % Selecting top %d percentile
				adj_tresh = prctile(Net_Adj(:), net_info.param_val);
			otherwise
				error('Unknown network type.');
		end
		Net_Adj = Net_Adj >= adj_tresh;
		Net_Adj(1:size(Net_Adj,1)+1:end) = 0; % Set diagonal to zero
		del_ind = sum(Net_Adj,1)==0; % Remove genes with no interactions
		Net_Adj(del_ind, :) = [];
		Net_Adj(:, del_ind) = [];
		Gene_Name(del_ind) = [];
		net_info.Net_Adj = Net_Adj;
		net_info.Gene_Name = Gene_Name;
end
end

function Dataset = getDataset(data_info, net_info, Valid_Gene_List)
	GeneInd = List2Index(Valid_Gene_List, data_info.Gene_Name);
	zData = zscore(data_info.Gene_Expression);
	Dataset.Gene_Expression = zData(data_info.CVInd, GeneInd);
	Dataset.Patient_Label = data_info.Patient_Label(data_info.CVInd);
	Dataset.Gene_Name = data_info.Gene_Name(GeneInd);
	if ~isequal(Dataset.Gene_Name, Valid_Gene_List), error('Gene name is inconsistent.'); end
	n_gene = size(Dataset.Gene_Expression, 2);
	
	if isempty(net_info), return; end
	switch net_info.net_name
		case 'Random'
			Net_Adj = rand(n_gene) > 0.98;
		case 'Corr'
			Net_Adj = double(corr(Dataset.Gene_Expression, 'Type', 'Spearman') > 0.6);
		case 'AbsCorr'
			Net_Adj = double(abs(corr(Dataset.Gene_Expression, 'Type', 'Spearman')) > 0.6);
		case {'DSN' 'CPN'}
			GeneInd = List2Index(Valid_Gene_List, net_info.Gene_Name);
			Net_Adj = net_info.Net_Adj(GeneInd, GeneInd);
		case {'BioGRID' 'CPDB' 'FunCoup-Complex' 'FunCoup-DOM' 'FunCoup-GIN' 'FunCoup-Metabolic' 'FunCoup-MEX' 'FunCoup-MIR' 'FunCoup-PEX' 'FunCoup-PFC' 'FunCoup-PHP' 'FunCoup-PPI' 'FunCoup-SCL' 'FunCoup-Signaling' 'FunCoup-TFB' 'HiC-ES' 'HiC-Full' 'HIN' 'HINT' 'HPRD' 'I2D' 'INstruct' 'iRefIndex' 'KEGG' 'MSigDB' 'Multinet' 'Reactome' 'String-CoEx' 'String-CoOc' 'String-Data' 'String-Exp' 'String-Fus' 'String-Neig' 'String-Txt' 'String' 'UniHI'}
			GMap = containers.Map(Dataset.Gene_Name, 1:n_gene);
			Net_Adj = zeros(n_gene);
			n_int = numel(net_info.net_cell{1});
			for ii=1:n_int
				showprogress(ii, n_int);
				if GMap.isKey(net_info.net_cell{1}{ii}) && GMap.isKey(net_info.net_cell{2}{ii})
					gi = GMap(net_info.net_cell{1}{ii});
					gj = GMap(net_info.net_cell{2}{ii});
					Net_Adj(gi, gj) = 1;
					Net_Adj(gj, gi) = 1;
				end
			end
		otherwise
			error('Unknown network.');
	end
	Net_Adj = double(max(Net_Adj, Net_Adj')>=1);
	Net_Adj(1:n_gene+1:end) = 1;
	if ~issymmetric(Net_Adj), error('Adj Matrix is not symetric.\n'); end
	Dataset.Net_Adj = Net_Adj;
end

function Ind_List = List2Index(Target_List, Population_List)
n_pop = numel(Population_List);
n_tar = numel(Target_List);
GMap = containers.Map(Population_List, 1:n_pop);

Ind_List = zeros(n_tar, 1);
for ti=1:n_tar
	Ind_List(ti) = GMap(Target_List{ti});
end
if ~isequal(Target_List, Population_List(Ind_List)), error(); end;
end