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
mi = find(ismember(clc_data.method_lst, {'TNMC','TNMCAd','TLEx','Lasso','GLasso'})); % 'GLassoAS','SGLasso','CFGLasso'
ni = find(ismember(clc_data.net_lst, {'AvgSynCrr-T10000','AbsCorr-T10000','STRING-T10000','KEGG-T10000','Random-T10000'}));
Method_lst = clc_data.method_lst(mi);
Net_lst = clc_data.net_lst(ni);
auc_mat = clc_data.Result_AUC(mi,ni,:,:);
[n_met, n_net, n_std, n_rep] = size(auc_mat);

%% Load baseline
load('../01_Pairwise_Evaluation_of_Genes/Baseline_AUCs/BA_CV01_TAgNMC_Random-T00020.mat', 'BaseLine_AUC');

%% Plotting
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(n_net);%[AdvancedColormap('cb', n_net/2); AdvancedColormap('mr', n_net/2)];
X_lbl = {};
step = 1;
% auc_mat = median(auc_mat, 4); % , 'omitnan'
auc_mat = mean(auc_mat, 3); % , 'omitnan'
for mi=1:n_met
	for ni=1:n_net
		auc_set = squeeze(auc_mat(mi, ni, :, :));% - BaseLine_AUC; % - squeeze(auc_mat(1, ni, :, :)); % - BaseLine_AUC(1:12); % 
		avg_set = mean(auc_set(:));
		std_set =  std(auc_set(:));
		if any(isnan(auc_set(:))), error(); end
		
		bar(step, avg_set, 'FaceColor', clr_map(ni,:));
		errorbarEx(step, avg_set, std_set, std_set, 2, 0.2, [0 0 0]);
		X_lbl{step, 1} = sprintf('%s', Method_lst{mi});
		step = step + 1;
	end
% 	text(step-n_met/2-1, 0.68, Net_lst{ni}, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top', ...
% 		'FontWeight', 'Bold', 'FontSize', 14, 'Rotation', 45);
	step = step + 1;
end
xlim([0 step-1]);
ylim([0.6 0.68]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(x) sprintf('%0.0f%%', x*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:step-2, 'XTickLabel', X_lbl, 'XTickLabelRotation', 20, ...
	'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10);
colormap(clr_map);
ylabel('Baseline Improvement', 'FontWeight', 'Bold');
legend(Net_lst, 'FontWeight', 'Bold', 'Location', 'NorthWest');

%% Saving
% output_name = sprintf('./Plots/S10_NOPs_Performance.pdf');
% set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [30 10], 'PaperPosition', [0 0 30 10]);
% print('-dpdf', '-r300', output_name);


