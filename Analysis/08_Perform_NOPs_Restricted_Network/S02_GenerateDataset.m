function ds_id = S02_GenerateDataset(cv_id, net_name)
%% ####
if ismac && ~exist('cv_id', 'var')
	fprintf('*** Warning!: Running on debug mode.\n');
	cv_id = 'Si01-Ri001';
	net_name = 'SyNet-T0500';
end

%% Initialization
cv_path = './CV_Files/';
dataset_path = './Dataset_Files/';
[~,~] = mkdir(dataset_path);
SEED_INFO=rng; Run_ID=double(SEED_INFO.Seed);
%rng shuffle % creates a different seed each time, but its not nessasary as seed is set in the job submission

%% Set Data ID
ds_lst = dir(sprintf([dataset_path '*_%s_%s_*'], cv_id, net_name));
if numel(ds_lst)>0
	ds_id = sprintf('%s_%s', cv_id, net_name);
	fprintf('[i] Dataset is already built for [%s]. Using [%s] ...\n', ds_id, ds_lst(1).name);
	return;
else
	ds_id = sprintf('%s_%s_%d', cv_id, net_name, Run_ID);
	dataset_name = sprintf([dataset_path 'DID_%s.mat'], ds_id);
	fprintf('Generating a dataset to be stored in [%s].\n', dataset_name);
end

%% Load CV info
cv_list = dir([cv_path '*' cv_id '*.mat']);
if numel(cv_list)~=1, error('Missing or duplicated CV info found. [%s]', strjoin({cv_list.name}, ', ')); end
cv_name = [cv_path cv_list.name];
fprintf('Loading CV info [%s] ...\n', cv_name);
load(cv_name, 'tr_info', 'te_info');
if isequal(tr_info.GEPath, te_info.GEPath) && any(tr_info.CVInd & te_info.CVInd)
	error('Some samples test set are used in training.');
end

%% Load Gene Expression
fprintf('Loading train expression data from [%s] ...\n', tr_info.GEPath);
data = load(tr_info.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
if ~isequal(data.Gene_Name, tr_info.Gene_Name), error(); end
tr_info.Gene_Expression = data.Gene_Expression;
tr_info.Patient_Label = (data.Patient_Label==1)*2-1;

fprintf('Loading test expression data from [%s] ...\n', te_info.GEPath);
data = load(te_info.GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
if ~isequal(data.Gene_Name, te_info.Gene_Name), error(); end
te_info.Gene_Expression = data.Gene_Expression;
te_info.Patient_Label = (data.Patient_Label==1)*2-1;
clear data
Valid_Gene_List = intersect(tr_info.Gene_Name, te_info.Gene_Name);

%% Load Network Info
fprintf('Loading [%s] network ...\n', net_name);
net_info = LoadNetworkInfo(net_name, tr_info, te_info);
Valid_Gene_List = intersect(Valid_Gene_List, net_info.Gene_Name);

%% Generate dataset
fprintf('Generating training set from [%s] ...\n', tr_info.GEName);
DatasetTr = UnifyData(tr_info, net_info, Valid_Gene_List);
DatasetTr.iCvPar = tr_info.iCvPar;
DatasetTr.GEName = tr_info.GEName;
DatasetTr.UsedSample = tr_info.CVInd;
DatasetTr.UsedStudy = tr_info.Study_Ind;
DatasetTr.GEPath = tr_info.GEPath;
DatasetTr.NetName = net_name;
DatasetTr.Net_info = rmfield(net_info, {'Net_Adj', 'Gene_Name'});

fprintf('Generating test set from [%s] ...\n', te_info.GEName);
DatasetTe = UnifyData(te_info, [], Valid_Gene_List);
DatasetTe.GEName = te_info.GEName;
DatasetTe.UsedSample = te_info.CVInd;
DatasetTe.UsedStudy = te_info.Study_Ind;
DatasetTe.GEPath = te_info.GEPath;

%% Save dataset
fprintf('Saving dataset in [%s].\n', dataset_name);
save(dataset_name, 'DatasetTr', 'DatasetTe');
end

function net_info = LoadNetworkInfo(net_name, tr_info, te_info)
dsn_path = '../01_Pairwise_Evaluation_of_Genes/Network_Files/';
name_part = regexp(net_name, '-', 'split');
net_info.net_name = name_part{1};
net_info.full_param = name_part{2};
net_info.param_type = name_part{2}(1);
net_info.param_val = str2double(name_part{2}(2:end));
switch net_info.net_name
	case 'Random'
		Gene_Name = tr_info.Gene_Name;
		Net_Adj = rand(numel(Gene_Name));
	case 'Corr'
		Net_Adj = corr(tr_info.Gene_Expression, 'Type', 'Spearman');
		Gene_Name = tr_info.Gene_Name;
	case 'AbsCorr'
		Net_Adj = abs(corr(tr_info.Gene_Expression, 'Type', 'Spearman'));
		Gene_Name = tr_info.Gene_Name;
	case 'SyNet'
		net_info.net_path = sprintf([dsn_path 'DSN_%sS%02d.mat'], net_info.net_name, te_info.Study_Ind);
		load(net_info.net_path, 'Net_Adj', 'Gene_Name');
	case {'CrMinSyn' 'CrSyn'}
		net_info.net_path = sprintf([dsn_path 'DSN_SyNetS%02d.mat'], te_info.Study_Ind);
		load(net_info.net_path, 'Pair_AUC', 'Gene_Name');
		n_gene = size(Pair_AUC,1);
		ind_auc = Pair_AUC(1:n_gene+1:end)';
		switch net_info.net_name
			case 'CrMinSyn'
				x_axis = bsxfun(@min, ind_auc, ind_auc');
			case 'CrSyn'
				x_axis = bsxfun(@(x,y) (x+y)/2, ind_auc, ind_auc');
		end
		ox = x_axis-min(x_axis(:));
		ox = ox./max(ox(:));
		
		pair_max = bsxfun(@max, ind_auc, ind_auc');
		y_axis = Pair_AUC./pair_max;
		oy = y_axis-min(y_axis(:));
		oy = oy./max(oy(:));
		
		z_data = load(tr_info.GEPath, 'Gene_Expression');
		z_axis = abs(corr(z_data.Gene_Expression(tr_info.CVInd,:), 'Type', 'Spearman'));
		z_axis(1:size(z_axis,1)+1:end) = 0;
		oz = z_axis-min(z_axis(:));
		oz = oz./max(oz(:));
		
		Net_Adj = single(-sqrt((ox-1).^2 + (oy-1).^2 + (oz-1).^2));
	case {'KEGG'}
		net_info.net_path = getPath(net_info.net_name);
		GSet_lst = regexp(fileread(net_info.net_path), '\n', 'split')';
		if strcmp(GSet_lst{end},''), GSet_lst(end)=[]; end
		n_gset = numel(GSet_lst);
		fprintf('Loading [%d] gene sets from [%s] ...\n', n_gset, net_info.net_path);
		for si=1:n_gset
			GSet_lst{si} = regexp(GSet_lst{si}, '\t', 'split');
		end
		Gene_Name = unique([GSet_lst{:}])';
		n_gene = numel(Gene_Name);
		GMap = containers.Map(Gene_Name, 1:n_gene);
		Net_Adj = zeros(n_gene);
		for si=1:n_gset
			grp_size = numel(GSet_lst{si});
			g_ind = zeros(grp_size, 1);
			for gi=1:grp_size
				g_ind(gi) = GMap(GSet_lst{si}{gi});
			end
			Net_Adj(g_ind, g_ind) = rand(grp_size);
		end
	case {'STRING'}
		net_info.net_path = getPath(net_info.net_name);
		fid = fopen(net_info.net_path, 'r');
		Header_lst = regexp(fgetl(fid), '\t', 'split');
		if numel(Header_lst)==2
			fprintf('No weight exists. Random selection of [%d] links.\n', net_info.param_val);
			net_cell = textscan(fid, '%s%s', 'Delimiter', '\t', 'ReturnOnError', 0);
			if ~feof(fid), error(); end
			n_lnk = numel(net_cell{1});
			rind = randperm(n_lnk, min([n_lnk net_info.param_val]));
			net_cell = [net_cell{1}(rind) net_cell{2}(rind)];
		else
			fprintf('Selecting of [%d] links from top weighted interactions.\n', net_info.param_val);
			net_cell = textscan(fid, '%s%s%d', net_info.param_val, 'Delimiter', '\t', 'ReturnOnError', 0);
		end
		fclose(fid);
		Gene_Name = unique(vertcat(net_cell{1:2}));
		n_gene = numel(Gene_Name);
		GMap = containers.Map(Gene_Name, 1:n_gene);
		Net_Adj = zeros(n_gene);
		n_int = numel(net_cell{1});
		fprintf('Forming the Adj matrix with [%d] genes: ', n_gene);
		for ii=1:n_int
			showprogress(ii, n_int);
			if GMap.isKey(net_cell{1}{ii}) && GMap.isKey(net_cell{2}{ii})
				gi = GMap(net_cell{1}{ii});
				gj = GMap(net_cell{2}{ii});
				Net_Adj(gi, gj) = n_int - ii + 1;
				Net_Adj(gj, gi) = n_int - ii + 1;
			end
		end
	otherwise
		error('Unknown network.');
end
Net_Adj = double(max(Net_Adj, Net_Adj'));
if ~issymmetric(Net_Adj), error('Adj Matrix is not symetric.\n'); end

%% Top selection
fprintf('Selecting top %d interactions.\n', net_info.param_val);
Net_Adj = Net_Adj - min(Net_Adj(:)); % Set minimum value to zero
Net_Adj = Net_Adj / max(Net_Adj(:)); % Set maximum value to one
Net_Adj(1:size(Net_Adj,1)+1:end) = 0; % Set diagonal to zero
scr_val = sort(Net_Adj(:), 'descend');
adj_tresh = scr_val(net_info.param_val);
Net_Adj(Net_Adj < adj_tresh) = 0;
fprintf('[%d] links are left in the network.\n', numel(nonzeros(Net_Adj(:))));

%% Node filtering
del_ind = sum(Net_Adj,1)==0; % Remove genes with no interactions
Net_Adj(del_ind, :) = [];
Net_Adj(:, del_ind) = [];
Gene_Name(del_ind) = [];
fprintf('[%d] genes are removed due to having no interactions.\n', sum(del_ind));

%% Storing
net_info.Net_Adj = Net_Adj;
net_info.Net_Threshold = adj_tresh;
net_info.Gene_Name = Gene_Name;
fprintf('[%d] genes are left in the network.\n', numel(Gene_Name));
end

function Dataset = UnifyData(data_info, net_info, Valid_Gene_List)
GeneInd = List2Index(Valid_Gene_List, data_info.Gene_Name);
zData = zscore(data_info.Gene_Expression);
Dataset.Gene_Expression = zData(data_info.CVInd, GeneInd);
Dataset.Patient_Label = data_info.Patient_Label(data_info.CVInd);
Dataset.Gene_Name = data_info.Gene_Name(GeneInd);
if ~isequal(Dataset.Gene_Name, Valid_Gene_List), error('Gene names are inconsistent.'); end
[n_sample, n_gene] = size(Dataset.Gene_Expression);
fprintf('After unification #Samples=%d, #Genes=%d are left.\n', n_sample, n_gene);

%% Network filtering
if isempty(net_info), return; end
GeneInd = List2Index(Valid_Gene_List, net_info.Gene_Name);
Net_Adj = net_info.Net_Adj(GeneInd, GeneInd);
% Net_Adj = double(max(Net_Adj, Net_Adj'));
Net_Adj(1:n_gene+1:end) = 0;
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
if ~isequal(Target_List, Population_List(Ind_List)), error(); end
end