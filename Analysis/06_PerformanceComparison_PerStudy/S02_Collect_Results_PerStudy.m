clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../_Utilities/');
res_path = './S01_Performance_PerStudy/';
Source_lst = {'ACES' 'METABRIC' 'TCGA'};
n_src = numel(Source_lst);

%% Load Study list
ge_name = getPath('SyNet');
fprintf('Loading expression from [%s] ...\n', ge_name);
load(ge_name, 'Patient_Info', 'Study_Index', 'Study_Name');
n_std = numel(Study_Name);

%% Main loop
auc_cell = cell(n_src, n_std);
for si=1:n_src
	for sj=1:n_std
		res_pattern = sprintf([res_path 'Res_LassoEx_%s_%s_SR0.70_*.mat'], Source_lst{si}, Study_Name{sj});
		res_lst = dir(res_pattern);
		n_res = numel(res_lst);
		fprintf('Found [%d] files for [%s] ...\n', n_res, res_pattern);
		for ri=1:n_res
			res_name = [res_path res_lst(ri).name];
			%fprintf('Loading [%s] ...\n', res_name);
			res_data = load(res_name, 'te_auc');
			auc_cell{si,sj} = [auc_cell{si,sj} res_data.te_auc];
		end
	end
end

%% Saving
save('Collected_Results.mat', 'Study_Name', 'Source_lst', 'auc_cell');
%{
%% Plotting
auc_mat = cellfun(@median, auc_cell);
figure('Position', [100 100 1400 800]);
img = imagesc(auc_mat);
set(img, 'AlphaData', ~isnan(auc_mat));
colormap(copper(10));
colorbar();
set(gca, 'XTick', 1:n_std, 'XTickLabel', Study_Name, 'XTickLabelRotation', 40, ...
		 'YTick', 1:n_src, 'YTickLabel', Source_lst, ...
		 'FontSize', 12, 'FontWeight', 'Bold');

%}