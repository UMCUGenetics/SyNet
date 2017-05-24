clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../_Utilities/');
res_path = './Results/';
Method_lst = {
	'NMC', 'LEx', ...
	'NMC-DSN-00100', 'NMC-DSN-00500', 'NMC-DSN-01000', 'NMC-DSN-05000', 'NMC-DSN-10000', 'NMC-DSN-20000', ...
	'LEx-DSN-00100', 'LEx-DSN-00500', 'LEx-DSN-01000', 'LEx-DSN-05000', 'LEx-DSN-10000', 'LEx-DSN-20000', ...
	'NMC-DsnGrp-00100', 'NMC-DsnGrp-00500', 'NMC-DsnGrp-01000', 'NMC-DsnGrp-05000', 'NMC-DsnGrp-10000', 'NMC-DsnGrp-20000', ...
	'LEx-DsnGrp-00100', 'LEx-DsnGrp-00500', 'LEx-DsnGrp-01000', 'LEx-DsnGrp-05000', 'LEx-DsnGrp-10000', 'LEx-DsnGrp-20000', ...
	'LEx-DsnIReg-00100', 'LEx-DsnIReg-00500', 'LEx-DsnIReg-01000', 'LEx-DsnIReg-05000', 'LEx-DsnIReg-10000', 'LEx-DsnIReg-20000', ...
	'LEx-FoldReg'
	};
n_met = numel(Method_lst);
n_rep = 5;

%% Load Study list
ge_name = getPath('SyNet');
fprintf('Loading studies from [%s] ...\n', ge_name);
load(ge_name, 'Study_Name');
n_study = 12; %numel(Study_Name);

%% Main loop
auc_mat = nan(n_met, n_study, n_rep);
for si=1:n_study
	for ri=1:n_rep
		for mi=1:n_met
			res_pattern = sprintf([res_path 'Res_%s_Si%02d_Ri%03d_*.mat'], Method_lst{mi}, si, ri);
			res_lst = dir(res_pattern);
			n_res = numel(res_lst);
			if n_res==0
				fprintf('[i] no result found for [%60s] ...\n', res_pattern);
				continue;
			end
			res_name = [res_path res_lst(1).name];
			fprintf('[%d] results: Reading [%45s] ...\n', n_res, res_lst(1).name);
			res_data = load(res_name, 'te_auc');
			auc_mat(mi,si,ri) = res_data.te_auc;
		end
	end
end
if any(isnan(auc_mat(:)))
	fprintf('Warning: Some nan values are detected.\n');
end

%% Saving
save('Collected_Results.mat', 'auc_mat', 'Method_lst', 'Study_Name');

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