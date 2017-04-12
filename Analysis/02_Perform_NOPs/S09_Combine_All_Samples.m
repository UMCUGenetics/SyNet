clc;
clear

%% Loading Top Pairs
data_lst = {'ACES' 'META' 'HAIBE'};
n_data = numel(data_lst);
for di=1:n_data
	%% Load gene expression
	fprintf('Loading [%s]\n', data_lst{di});
	GeneExpression_Path = getPath(data_lst{di});
	data = load(GeneExpression_Path, 'Gene_Expression', 'Patient_Label', 'Study_Index', 'Study_Name', 'Gene_Name');
	cls_freq = histcounts(data.Patient_Label, 0:2);
	fprintf('Loaded [%d/%d] (%0.1f, %0.1f)%% classes of patients.\n', histcounts(data.Patient_Label), cls_freq*100/sum(cls_freq));
	if ~isfield(data, 'Study_Name')
		data.Study_Index = ones(size(data.Patient_Label));
		data.Study_Name{1} = data_lst{di};
	end
	CMB_data(di) = data;
	clear data
end

%% Unify gene sets
Gene_Name = CMB_data(1).Gene_Name;
for di=2:n_data
	Gene_Name = intersect(Gene_Name, CMB_data(di).Gene_Name);
end

%% Select Anchor gene expressions
Gene_Expression = [];
Patient_Label = [];
Study_Index = [];
Study_Name = [];
for di=1:n_data
	ind = List2Index(Gene_Name, CMB_data(di).Gene_Name);
	offset = numel(Study_Name);
	Gene_Expression = [Gene_Expression; zscore(CMB_data(di).Gene_Expression(:, ind))];
	Patient_Label = [Patient_Label; CMB_data(di).Patient_Label];
	Study_Index = [Study_Index; offset + CMB_data(di).Study_Index];
	Study_Name = [Study_Name; strcat(data_lst{di},'-', CMB_data(di).Study_Name)];
end
n_sample = size(Gene_Expression, 1);

%% Correct the index for missing studies
[std_uid, ~, std_nid] = unique(Study_Index, 'stable');
Study_Name = Study_Name(std_uid);
Study_Index = std_nid;

%% Compute sample correlation
sam_cr = corr(Gene_Expression', 'Type', 'Spearman');
sam_cr(1:n_sample+1:end) = 0;
high_cr = find(max(abs(sam_cr))>0.8);
sel_cr = sam_cr(high_cr, high_cr);
study_lbl = Study_Index(high_cr);
sel_cr(:, end+1) = study_lbl/max(study_lbl);
sel_cr(end+1,1:end-1) = sel_cr(:, end);
Study_Name(Study_Index(high_cr))
if ismac
	img_h = imagesc(sel_cr);
	set(img_h, 'AlphaData', abs(sel_cr)>=0.5);
	colormap(flipud(jet(30)));
end

%% Save gene expression
out.Gene_Expression = Gene_Expression;
out.Patient_Label = Patient_Label;
out.Study_Index = Study_Index;
out.Study_Name = Study_Name;
out.Gene_Name = Gene_Name;
save('./CMB_ACES-META-HAIBE.mat', '-struct', 'out');
fprintf('Combination finished.\n');

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