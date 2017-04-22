clc;
clear;

%% Initialization
n_study = 15;
addpath('../../../../Useful_Sample_Codes/getTop');
study_name = cellstr(num2str((1:n_study)','./Network_Files/DSN_SyNetS%02d.mat'));
study_name{n_study} = './Network_Files/DSN_SyNet.mat';
n_top = 10000;

%% Main loop
Top_Gene = cell(n_study, 1);
for si=1:n_study
	fprintf('Loading [%s]...\n', study_name{si});
	load(study_name{si}, 'Net_Adj');
	n_gene = size(Net_Adj,1);
	fprintf('Sorting ...\n');
	[scr_val, scr_ind] = sort(Net_Adj(:), 'Descend');
	%[i,j] = find(Net_Adj>=scr_val(n_top));
	%[Top_Ind, Top_Freq] = getTop([i;j], 100);
	[i, j] = ind2sub([n_gene, n_gene], scr_ind(1:n_top));
	Top_Ind = unique([i j]', 'Stable'); Top_Ind = Top_Ind(1:500);
	Top_Gene{si,1} = Top_Ind;
end

%% Measure overlap
ovl_mat = zeros(n_study);
for si=1:n_study
	for sj=1:n_study
		ovl_mat(si,sj) = numel(intersect(Top_Gene{si}, Top_Gene{sj}));
	end
end

%% Plot overlap
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
Study_Name{end+1} = 'All';
imagesc(ovl_mat);
set(gca, 'XTick', 1:n_study, 'XTickLabel', Study_Name, 'XTickLabelRotation', 25, ...
		 'YTick', 1:n_study, 'YTickLabel', Study_Name, 'YTickLabelRotation', 25, 'FontWeight', 'Bold', 'FontSize', 14);
colormap(copper(10));
clr_h = colorbar();
ylabel(clr_h, 'Number of overlapping genes');
	 