clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../_Utilities/');
res_path = './Results_Files/';
Method_lst = {'TLEx' 'LExAG' 'Lasso' 'iPark' 'RI-iPark' 'Feral' 'RI-Feral' 'AS-Feral'};
n_met = numel(Method_lst);
n_pair = 10000;
Net_lst = strcat({'Random-T' 'Corr-T' 'AbsCorr-T' 'SyNet-T'}, sprintf('%05d', n_pair));
n_net = numel(Net_lst);
n_study = 12; %numel(Study_Name);
n_rep = 10;
n_MAX_feat = 500;

%% Load Study list
ge_name = getPath('SyNet');
fprintf('Loading studies from [%s] ...\n', ge_name);
load(ge_name, 'Study_Name');

%% Main loop
auc_mat = nan(n_met, n_net, n_study, n_rep);
for si=1:n_study
	fprintf('Study %d ...\n', si);
	for ri=1:n_rep
		for mi=1:n_met
			for ni=1:n_net
				res_pattern = sprintf([res_path 'DID_Si%02d-Ri%03d_%s_*_MSN-%03d_MTN-%s.mat'], si, ri, Net_lst{ni}, n_MAX_feat, Method_lst{mi});
				res_lst = dir(res_pattern);
				n_res = numel(res_lst);
				if n_res==0
					fprintf('********* no result found for [%60s] ...\n', res_pattern);
					continue;
				end
				res_name = [res_path res_lst(1).name];
				%fprintf('[%d] results: Reading [%45s] ...\n', n_res, res_lst(1).name);
				res_data = load(res_name, 'te_auc');
				auc_mat(mi,ni,si,ri) = res_data.te_auc;
			end
		end
	end
end
if any(isnan(auc_mat(:)))
	fprintf('Warning: Some nan values are detected.\n');
end

%% Saving
sav_name = sprintf('Collected_Results_MSN-%03d.mat', n_MAX_feat);
save(sav_name, 'auc_mat', 'Method_lst', 'Study_Name', 'Net_lst');

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