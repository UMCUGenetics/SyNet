clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/getRankAUC/');
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/FisherExactTest/');
gs_path = './Enrichment_Sets/';
ntp_path = './Net_TopPairs/';
[~,~] = mkdir(ntp_path);
Net_lst = {'TTest', 'Synergy', 'Corr', 'SyNet', 'Random'}; %  , 'STRING', 'KEGG'
n_net = numel(Net_lst);
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name');
n_gene = numel(Gene_Name);
n_top = n_gene;

%% Load Top genes in correlation
Mrk_lst = cell(n_net, 1);
for ni=1:n_net
	sav_name = sprintf([ntp_path 'TopPair_%s_%05d.mat'], Net_lst{ni}, n_top);
	if ~exist(sav_name, 'File')
		switch Net_lst{ni}
			case 'TTest'
				Top_Gene = getTopTTest(n_top);
			case 'Corr'
				Top_Gene = getTopCorr(n_top);
			case 'Synergy'
				Top_Gene = getTopSynergy(n_top);
			case 'SyNet'
				Top_Gene = getTopSyNet(n_top);
			case 'STRING'
				Top_Gene = getTopNet('../../Networks/STRING/STRING_combined_score_GN_WithWeights.txt', n_top);
			case 'KEGG'
				Top_Gene = getTopSet('../../Networks/KEGG/KEGG.txt', n_top);
			case 'Random'
				rind = randperm(numel(Gene_Name), n_top);
				Top_Gene = Gene_Name(rind);
		end
		save(sav_name, 'Top_Gene');
	else
		fprintf('Loading [%s] ... \n', sav_name);
		load(sav_name);
	end
	Mrk_lst{ni} = Top_Gene;
end

%% Load gene set
gset_file = dir([gs_path 'RefSet_*.csv']);
n_gset = numel(gset_file);
gset_list = cell(n_gset,1);
gset_name = cell(n_gset,1);
for si=1:n_gset
    gset_list{si} = regexp(fileread([gs_path gset_file(si).name]), '\n', 'split')';
	gset_list{si} = intersect(gset_list{si}, Gene_Name);
    file_info = regexp(gset_file(si).name, '[_\.]', 'split');
    gset_name{si} = file_info{2};
end

%% Measure enrichment score
auc_mat = zeros(n_gset, n_net);
for ni=1:n_net
	Mrk_set = Mrk_lst{ni};
	%rest_set = setdiff(Gene_Name, Mrk_set);
	%n_rest = numel(rest_set);
	%GSet = Mrk_set; %rest_set(randperm(n_rest))];
	for si=1:n_gset
		auc_mat(si, ni) = getRankAUC(gset_list{si}, Mrk_set);
	end
end

%% Plotting
close all
figure('Position', [100 100 1600 500]);
bar(auc_mat, 'grouped');
legend(Net_lst, 'Location', 'NorthEast', 'Orientation', 'Vertical');

yTick = 0:0.05:1;
yTickLabel = arrayfun(@(x) num2str(x), yTick, 'UniformOutput', false);
set(gca, 'XTick', 1:n_gset, 'XTicklabel', gset_name, 'XTicklabelRotation', 20, 'yTick', yTick, 'yticklabel', yTickLabel, ...
	'ygrid', 'on', 'fontweight', 'bold', 'fontsize', 16);
xlim([0 n_gset+1]);
ylim([0.5 0.7]);
clr_map = jet(n_net);
% clr_map(2,:) = [];
colormap(clr_map);

%% Save plot
% print('-dpdf', '-r300', ['./Plots/Gene_Enrichment.pdf']);

%%//////////////////////////////////////////////////////////////////////////////////
function Top_Pair = getTopPairs(Net_Adj)
n_gene = size(Net_Adj, 1);
n_total = n_gene*(n_gene-1)/2;
Pair_Info = zeros(n_total, 3, 'single');
[Pair_Info(:,1), Pair_Info(:,2)] = find(triu(ones(n_gene), 1));
pair_ind = sub2ind([n_gene n_gene], Pair_Info(:,1), Pair_Info(:,2));
Pair_Info(:,3) = Net_Adj(pair_ind);
[~, sid] = sort(Pair_Info(:,3), 'Descend');
Top_Pair = Pair_Info(sid, :);
end

function Top_Gene = getTopTTest(n_top)
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Expression', 'Patient_Label', 'Gene_Name');
zData = zscore(Gene_Expression);
Patient_Label = (Patient_Label==1)*2-1;
n_gene = numel(Gene_Name);
pv_vec = ttest2Ex(zData, Patient_Label);
[top_val, top_ind] = sort(-log10(pv_vec), 'Descend');
Top_Gene = Gene_Name(top_ind(1:n_top));
end

function Top_Gene = getTopCorr(n_top)
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Expression', 'Gene_Name');
Net_Adj = abs(corr(Gene_Expression, 'Type', 'Spearman'));
Top_Pair = getTopPairs(Net_Adj);
top_ind = unique(Top_Pair(:,1:2)', 'Stable');
Top_Gene = Gene_Name(top_ind(1:n_top));
end

function Top_Gene = getTopSyNet(n_top)
net_path = '../01_Pairwise_Evaluation_of_Genes/Network_Files/DSN_SyNet.mat';
ge_path = '../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat';
load(net_path, 'Pair_AUC', 'Gene_Name');
n_gene = size(Pair_AUC,1);
ind_auc = Pair_AUC(1:n_gene+1:end)';
x_axis = bsxfun(@(x,y) (x+y)/2, ind_auc, ind_auc');
ox = x_axis-min(x_axis(:));
ox = ox./max(ox(:));

pair_max = bsxfun(@max, ind_auc, ind_auc');
y_axis = Pair_AUC./pair_max;
oy = y_axis-min(y_axis(:));
oy = oy./max(oy(:));

GE_Data = load(ge_path, 'Gene_Expression', 'Gene_Name');
if ~isequal(Gene_Name, GE_Data.Gene_Name), error(); end
z_axis = abs(corr(GE_Data.Gene_Expression, 'Type', 'Spearman'));
z_axis(1:size(z_axis,1)+1:end) = 0;
oz = z_axis-min(z_axis(:));
oz = oz./max(oz(:));

Net_Adj = single(-sqrt((ox-1).^2 + (oy-1).^2 + (oz-1).^2));
Top_Pair = getTopPairs(Net_Adj);
top_ind = unique(Top_Pair(:,1:2)', 'Stable');
Top_Gene = Gene_Name(top_ind(1:n_top));
end

function Top_Gene = getTopSynergy(n_top)
net_path = '../01_Pairwise_Evaluation_of_Genes/Network_Files/DSN_SyNet.mat';
ge_path = '../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat';
load(net_path, 'Pair_AUC', 'Gene_Name');
n_gene = size(Pair_AUC,1);

ind_auc = Pair_AUC(1:n_gene+1:end)';
pair_max = bsxfun(@max, ind_auc, ind_auc');
y_axis = Pair_AUC./pair_max;
oy = y_axis-min(y_axis(:));
oy = oy./max(oy(:));

GE_Data = load(ge_path, 'Gene_Expression', 'Gene_Name');
if ~isequal(Gene_Name, GE_Data.Gene_Name), error(); end
z_axis = abs(corr(GE_Data.Gene_Expression, 'Type', 'Spearman'));
z_axis(1:size(z_axis,1)+1:end) = 0;
oz = z_axis-min(z_axis(:));
oz = oz./max(oz(:));

Net_Adj = single(-sqrt((oy-1).^2 + (oz-1).^2));
Top_Pair = getTopPairs(Net_Adj);
top_ind = unique(Top_Pair(:,1:2)', 'Stable');
Top_Gene = Gene_Name(top_ind(1:n_top));
end

function Top_Gene = getTopNet(net_path, n_top)
n_top_pairs = n_top*100;
fid = fopen(net_path, 'r');
Header_lst = regexp(fgetl(fid), '\t', 'split');
if numel(Header_lst)==2
	fprintf('No weight exists. Random selection of [%d] links.\n', n_top_pairs);
	net_cell = textscan(fid, '%s%s', 'Delimiter', '\t', 'ReturnOnError', 0);
	if ~feof(fid), error(); end
	n_lnk = numel(net_cell{1});
	rind = randperm(n_lnk, min([n_lnk net_info.param_val]));
	net_cell = [net_cell{1}(rind) net_cell{2}(rind)];
else
	fprintf('Selection of [%d] links from top weighted interactions.\n', n_top_pairs);
	net_cell = textscan(fid, '%s%s%d', n_top_pairs, 'Delimiter', '\t', 'ReturnOnError', 0);
end
fclose(fid);
Gene_Name = unique(vertcat(net_cell{1:2}));
n_gene = numel(Gene_Name);
GMap = containers.Map(Gene_Name, 1:n_gene);
Net_Adj = zeros(n_gene);
n_int = numel(net_cell{1});
fprintf('Forming the Adj matrix with [%d] genes\n', n_gene);
for ii=1:n_int
	if GMap.isKey(net_cell{1}{ii}) && GMap.isKey(net_cell{2}{ii})
		gi = GMap(net_cell{1}{ii});
		gj = GMap(net_cell{2}{ii});
		Net_Adj(gi, gj) = n_int - ii + 1;
		Net_Adj(gj, gi) = n_int - ii + 1;
	end
end

Top_Pair = getTopPairs(Net_Adj);
top_ind = unique(Top_Pair(:,1:2)', 'Stable');
Top_Gene = Gene_Name(top_ind(1:min([n_top numel(top_ind)])));
end

function Top_Gene = getTopSet(set_path, n_top)
GSet_lst = regexp(fileread(set_path), '\n', 'split')';
if strcmp(GSet_lst{end},''), GSet_lst(end)=[]; end
n_gset = numel(GSet_lst);
fprintf('Loading [%d] gene sets from [%s] ...\n', n_gset, set_path);
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

Top_Pair = getTopPairs(Net_Adj);
top_ind = unique(Top_Pair(:,1:2)', 'Stable');
Top_Gene = Gene_Name(top_ind(1:min([n_top numel(top_ind)])));
end