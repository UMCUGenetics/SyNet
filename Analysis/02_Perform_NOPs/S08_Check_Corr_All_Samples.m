clc;
clear

%% Loading Top Pairs
data_lst = {'ACES' 'META' 'HAIBE'};
Sample_Name = {};
for di=1:numel(data_lst)
	%% Load gene expression
	fprintf('Loading [%s]\n', data_lst{di});
	GeneExpression_Path = getPath(data_lst{di});
	out = load(GeneExpression_Path, 'Gene_Expression', 'Gene_Name');
	data(di).Gene_Expression = zscore(out.Gene_Expression);
	data(di).Gene_Name = out.Gene_Name;
	data(di).GE_Name = data_lst{di};
	Sample_Name = [Sample_Name; repmat(data_lst(di), size(out.Gene_Expression,1), 1)];
end

%% Unify gene sets
Gene_Name = intersect(data(1).Gene_Name, data(2).Gene_Name);
Gene_Name = intersect(Gene_Name, data(3).Gene_Name);

ind = List2Index(Gene_Name, data(1).Gene_Name);
Gene_Expression = data(1).Gene_Expression(:, ind);

ind = List2Index(Gene_Name, data(2).Gene_Name);
Gene_Expression = [Gene_Expression; data(2).Gene_Expression(:, ind)];

ind = List2Index(Gene_Name, data(3).Gene_Name);
Gene_Expression = [Gene_Expression; data(3).Gene_Expression(:, ind)];
[n_sample, n_gene] = size(Gene_Expression);

%% Compute sample correlation
sam_cr = corr(Gene_Expression', 'Type', 'Spearman');
sam_cr(1:n_sample+1:end) = 0;
high_cr = find(max(abs(sam_cr))>0.8);
sel_cr = sam_cr(high_cr, high_cr);
sel_cr(:, end+1) = ismember(Sample_Name(high_cr), data_lst{2})*1 + ismember(Sample_Name(high_cr), data_lst{3})*-1;
sel_cr(end+1,1:end-1) = sel_cr(:, end);
img_h = imagesc(sel_cr);
set(img_h, 'AlphaData', abs(sel_cr)>=0.6);
colormap(flipud(jet(15)));

%% Check correlating patients
max_cr = max(sam_cr, [], 2);
in = find(max_cr>0.8);
Sample_Name(in)

function path_str = getPath(data_name)
switch data_name
	case 'ACES'
		path_str = '../../../Dataset/ACES/GE.mat';
	case 'GBM'
		path_str = '../../../Dataset/Cancer_Genomics_Browser/CGB_GBM.mat';
	case 'GBMFl'
		path_str = '../95_Filter_GBM/CGB_GBM_FreqFiltered.mat';
	case 'META'
		path_str = '../94_Pairwise_Evaluation_of_Genes_METABRIC/METABRIC_Filtered.mat';
	case 'HAIBE'
		path_str = '../96_Pairwise_Evaluation_of_Genes_Haibe/Haibe_Filtered.mat';
		%path_str = '../96_Pairwise_Evaluation_of_Genes_Haibe/Haibe_Filtered_New.mat';
		%path_str = '../96_Pairwise_Evaluation_of_Genes_Haibe/Haibe_Filtered_Old.mat';
end
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