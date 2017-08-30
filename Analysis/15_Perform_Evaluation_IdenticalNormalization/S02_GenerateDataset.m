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
if ~ismember(net_info.param_type, {'P','G'}), error(); end
net_info.param_val = str2double(name_part{2}(2:end));
if strcmp(net_info.param_type, 'P')
	MAX_N_PAIR = net_info.param_val;
else
	MAX_N_PAIR = inf;
end
switch net_info.net_name
	case 'None'
        fprintf('No network is chosen. Ignoring network preparation...\n');
		net_info.Gene_Name = tr_info.Gene_Name;
		net_info.Net_Adj = zeros(numel(net_info.Gene_Name), 'single');
        net_info.net_path = '';
        return;
    case 'Random'
		Gene_Name = tr_info.Gene_Name;
		Net_Adj = rand(numel(Gene_Name));
	case 'Corr'
		Net_Adj = corr(tr_info.Gene_Expression, 'Type', 'Spearman');
		Gene_Name = tr_info.Gene_Name;
	case 'AbsCorr'
		Net_Adj = abs(corr(tr_info.Gene_Expression, 'Type', 'Spearman'));
		Gene_Name = tr_info.Gene_Name;
	case {'KEGG', 'MSigDB'}
		net_info.net_path = getPath(net_info.net_name);
		GSet_lst = regexp(fileread(net_info.net_path), '\n', 'split')';
		if strcmp(GSet_lst{end},''), GSet_lst(end)=[]; end
		n_gset = numel(GSet_lst);
		fprintf('Loading [%d] gene sets from [%s] ...\n', n_gset, net_info.net_path);
		for si=1:n_gset
			GSet_lst{si} = regexp(GSet_lst{si}, '\t', 'split');
		end
		Gene_Name = unique([GSet_lst{:}])';
		fprintf('[i] Network contains [%d] genes before filtering.\n', numel(Gene_Name));
		Gene_Name = intersect(intersect(Gene_Name, tr_info.Gene_Name), te_info.Gene_Name); %% Filter extra genes
		fprintf('[i] Network contains [%d] genes after filtering.\n', numel(Gene_Name));
		n_gene = numel(Gene_Name);
		GMap = containers.Map(Gene_Name, 1:n_gene);
		Net_Adj = zeros(n_gene, 'single');
		fprintf('Adding genes to network: ');
		for si=1:n_gset
			showprogress(si, n_gset, 20);
			grp_size = numel(GSet_lst{si});
			g_ind = zeros(grp_size, 1);
			for gi=1:grp_size
				if GMap.isKey(GSet_lst{si}{gi}) 
					g_ind(gi) = GMap(GSet_lst{si}{gi});
				end
			end
			g_ind = nonzeros(g_ind);
			Net_Adj(g_ind, g_ind) = rand(numel(g_ind));
		end
		clear GSet_lst
	case {'STRING','HPRD','I2D'}
		net_info.net_path = getPath(net_info.net_name);
		fid = fopen(net_info.net_path, 'r');
		Header_lst = regexp(fgetl(fid), '\t', 'split');
		if numel(Header_lst)==2
			fprintf('No weight exists. Random selection of [%d] links.\n', MAX_N_PAIR);
			net_cell = textscan(fid, '%s%s', 'Delimiter', '\t', 'ReturnOnError', 0);
			if ~feof(fid), error(); end
			n_lnk = numel(net_cell{1});
			rind = randperm(n_lnk, min([n_lnk MAX_N_PAIR]));
			net_cell = {net_cell{1}(rind) net_cell{2}(rind)};
		else
			fprintf('Selecting of [%d] links from top weighted interactions.\n', MAX_N_PAIR);
			net_cell = textscan(fid, '%s%s%d', MAX_N_PAIR, 'Delimiter', '\t', 'ReturnOnError', 0);
		end
		fclose(fid);
		Gene_Name = unique(vertcat(net_cell{1:2}));
		fprintf('[i] Network contains [%d] genes before filtering.\n', numel(Gene_Name));
		Gene_Name = intersect(intersect(Gene_Name, tr_info.Gene_Name), te_info.Gene_Name); %% Filter extra genes
		fprintf('[i] Network contains [%d] genes after filtering.\n', numel(Gene_Name));
		n_gene = numel(Gene_Name);
		GMap = containers.Map(Gene_Name, 1:n_gene);
		Net_Adj = zeros(n_gene, 'single');
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
		clear net_cell
	otherwise
		net_info.net_path = sprintf([dsn_path 'DSN_SyNetS%02d.mat'], te_info.Study_Ind);
		load(net_info.net_path, 'Pair_AUC', 'Gene_Name');
		n_gene = size(Pair_AUC,1);
		ind_auc = Pair_AUC(1:n_gene+1:end)';
		Pair_Dist = zeros(n_gene);
		for ni=1:3:numel(net_info.net_name)
			nn_part = net_info.net_name(ni:ni+2);
			fprintf('Computing [%s] ... \n', nn_part);
			switch nn_part
				case 'Avg'
					ax_avg = bsxfun(@(x,y) (x+y)/2, ind_auc, ind_auc');
					Pair_Dist = Pair_Dist + (oscore(ax_avg)-1).^2;
					clear ax_avg
				case 'ACr'
					ge_data = load(tr_info.GEPath, 'Gene_Expression');
					ax_crr = abs(corr(ge_data.Gene_Expression(tr_info.CVInd,:), 'Type', 'Spearman'));
					ax_crr(1:size(ax_crr,1)+1:end) = 0;
					Pair_Dist = Pair_Dist + (oscore(ax_crr)-1).^2;
					clear ge_data ax_crr
				case 'Syn'
					pair_max = bsxfun(@max, ind_auc, ind_auc');
					Pair_Dist = Pair_Dist + (oscore(Pair_AUC./pair_max)-1).^2;
					clear pair_max
				otherwise
					error('Unknown error.');
			end
		end
		Net_Adj = single(-sqrt(Pair_Dist));
		clear NetScr
end
clear tr_info te_info

fprintf('Making sure the network is symetric ...\n');
Net_Adj = max(Net_Adj, Net_Adj');
if ~issymmetric(Net_Adj), error('Adj Matrix is not symetric.\n'); end

%% Top selection
fprintf('Normalizing edge weights between [0 1] ...\n');
Net_Adj = Net_Adj - min(Net_Adj(:)); % Set minimum value to zero
Net_Adj = Net_Adj / max(Net_Adj(:)); % Set maximum value to one
Net_Adj(1:size(Net_Adj,1)+1:end) = 0; % Set diagonal to zero
Net_Eps = Net_Adj + rand(size(Net_Adj,1))*1e-10; % Add small variation to make sure links with same weight do not exists
[scr_val, scr_ind] = sort(Net_Eps(:), 'Descend');
if strcmp(net_info.param_type, 'P')
	clear scr_ind
	fprintf('Selecting top %d interactions.\n', MAX_N_PAIR);
	adj_tresh = scr_val(MAX_N_PAIR);
	Net_Adj(Net_Eps < adj_tresh) = 0;
	net_info.Net_Threshold = adj_tresh;
	clear Net_Eps scr_val
else
	clear Net_Eps scr_val
	n_gene = numel(Gene_Name);
	MAX_N_GENE = min([n_gene net_info.param_val]);
    fprintf('Selecting top %d genes.\n', MAX_N_GENE);
    top_ind = zeros(numel(scr_ind), 2, 'uint16');
    [top_ind(:,1), top_ind(:,2)] = ind2sub([n_gene n_gene], scr_ind);
    if any(top_ind(:)>=65535), error('Not implemented for such a large number of genes!'); end
    gin_map = containers.Map('KeyType', 'Double', 'ValueType', 'Double');
    for rest_ind=1:numel(scr_ind)
        if ~gin_map.isKey(top_ind(rest_ind,1)), gin_map(top_ind(rest_ind,1)) = 0; end
        if ~gin_map.isKey(top_ind(rest_ind,2)), gin_map(top_ind(rest_ind,2)) = 0; end
        if gin_map.Count >= MAX_N_GENE
            break;
        end
    end
    Net_Adj(scr_ind(rest_ind+1:end)) = 0;
    Net_Adj = max(Net_Adj, Net_Adj');
    clear top_ind scr_ind
end
fprintf('[%d] links are left in the network.\n', numel(nonzeros(Net_Adj(:))));

%% Node filtering
fprintf('Removing genes with no interactions ...\n');
del_ind = sum(Net_Adj,1)==0;
Net_Adj(del_ind, :) = [];
Net_Adj(:, del_ind) = [];
Gene_Name(del_ind) = [];
fprintf('[%d] genes are removed due to having no interactions.\n', sum(del_ind));

%% Storing
net_info.Net_Adj = double(Net_Adj);
net_info.Gene_Name = Gene_Name;
fprintf('[%d] genes and [%d] links are left in the network.\n', numel(Gene_Name), numel(nonzeros(Net_Adj(:))));
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

function lst = oscore(lst)
lst = lst - min(lst(:));
lst = lst ./ max(lst(:));
end