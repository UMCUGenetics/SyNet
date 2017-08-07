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
mi = [1 3 5 9 12];
ni = [7 12:16];
Method_lst = clc_data.method_lst;%(mi);
Net_lst = clc_data.net_lst;%(ni);
auc_mat = clc_data.Result_AUC;%(mi,ni,:,:);
n_met = numel(Method_lst);
n_net = numel(Net_lst);
n_study = 12;
n_rep = 10;

%% Collect AUC
auc_std_med = nan(n_met, n_net);
for mi=1:n_met
	for ni=1:n_net
		auc_std = std(squeeze(auc_mat(mi, ni, :, :)),0,2);
		auc_std_med(mi, ni) = median(auc_std, 'omitnan');
	end
end

%% Plot
img_h = imagesc(auc_std_med);
set(img_h, 'AlphaData', ~isnan(auc_std_med));
set(gca, 'YTick', 1:n_met, 'YTickLabel', Method_lst, 'XTickLabelRotation', 20, ...
		 'XTick', 1:n_net, 'XTickLabel', Net_lst, 'FontWeight', 'Bold', 'FontSize', 10);
colorbar();
colormap(jet(5));
caxis([0 0.03]);
%% Saving
output_name = sprintf('./Plots/S08_PerfStabilityMatrix_Per-Method-Net.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [15 10], 'PaperPosition', [0 0 15 10]);
print('-dpdf', '-r300', output_name);

	