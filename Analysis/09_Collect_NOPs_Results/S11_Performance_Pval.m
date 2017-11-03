clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
res_name = 'Collected_NOP_CV01_AUCs.mat';
fprintf('Loading results from [%s] ...\n', res_name);
clc_data = load(res_name);
mi = [3 5 12];
ni = [7 12]; %  14 16
Method_lst = clc_data.method_lst(mi);
Net_lst = clc_data.net_lst(ni);
auc_mat = clc_data.Result_AUC(mi, ni, :, :);
n_met = numel(Method_lst);
n_net = numel(Net_lst);
n_std = 12;
n_rep = 10;

%% Plotting
step = 1;
X_lbl = {};
auc_cell = {};
% auc_mat = mean(auc_mat, 4);
for ni=1:n_net
	for mi=1:n_met
		auc_cell{step,1} = squeeze(auc_mat(mi, ni, :));
		X_lbl{step, 1} = sprintf('%s,%s', Method_lst{mi}, Net_lst{ni});
		step = step + 1;
	end
end
n_elm = numel(auc_cell);

scr_mat = zeros(n_elm);
for ei=1:n_elm
	for ej=1:n_elm
		[h, pval, ci, stats] = ttest(auc_cell{ei}, auc_cell{ej});
		if stats.tstat>0
			scr_mat(ei,ej) = -log10(pval);
		else
			scr_mat(ei,ej) = log10(pval);
		end
	end
end

%% Plotting
close all
img_h = imagesc(scr_mat);
set(img_h, 'AlphaData', ~isnan(scr_mat));
clr_h = colorbar();
clr_map = AdvancedColormap('rwb', 15);
colormap(clr_map);
ylabel(clr_h, '-Log10(pval)');
set(gca, 'XTick', 1:n_elm, 'XTickLabel', X_lbl, 'XTickLabelRotation', 35, ...
		'YTick', 1:n_elm, 'YTickLabel', X_lbl, 'YDir', 'Normal', ...
	'FontWeight', 'Bold', 'FontSize', 7);
caxis([-8 8]);

%% Saving
output_name = sprintf('./Plots/S11_MetNet_AUCPval.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [30 25], 'PaperPosition', [0 0 30 25]);
print('-dpdf', '-r300', output_name);
fprintf('\n');

