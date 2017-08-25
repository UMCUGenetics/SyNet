clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
cv_ind = 1;
res_name = sprintf('Collected_NOP_CV%02d_Markers.mat', cv_ind);
fprintf('Loading results from [%s] ...\n', res_name);
clc_data = load(res_name);
if any(isnan(clc_data.Result_AUC(:)))
	fprintf('Warning, some NaN value were detected.\n'); 
end

%% Select methods and networks
mi = getInd(clc_data.method_lst, {'TNMC','TNMCAd','TLEx','Lasso','GLasso'}); % 'GLassoAS','SGLasso','GLasso','CFGLasso'
ni = getInd(clc_data.net_lst, {
	'AbsCorr-G00500','AbsCorr-P10000','AvgSyn-G00500','AvgSyn-P10000','AvgSynACr-G00500','AvgSynACr-P10000', ...
	'HPRD-G00500','HPRD-P10000','I2D-G00500','I2D-P10000','KEGG-G00500','KEGG-P10000', ...
	'MSigDB-G00500','MSigDB-P10000','Random-G00500','Random-P10000','STRING-G00500','STRING-P10000' ...
});
% mi=1:numel(clc_data.method_lst);
% ni=1:numel(clc_data.net_lst);

Method_lst = clc_data.method_lst(mi);
Net_lst = clc_data.net_lst(ni);
Result_AUC = clc_data.Result_AUC(mi,ni,:,:);
[n_met, n_net, n_std, n_rep] = size(Result_AUC);

%% Plotting
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(n_net);%[AdvancedColormap('cb', n_net/2); AdvancedColormap('mr', n_net/2)];
X_lbl = {};
step = 1;
% Result_AUC = median(Result_AUC, 3); % , 'omitnan'
Result_AUC = mean(Result_AUC, 3); % , 'omitnan'
for mi=1:n_met
	for ni=1:n_net
		auc_set = squeeze(Result_AUC(mi, ni, :, :));% - BaseLine_AUC; % - squeeze(auc_mat(1, ni, :, :)); % - BaseLine_AUC(1:12); % 
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
ylim([58 70]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(x) sprintf('%0.0f%%', x), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:step-2, 'XTickLabel', X_lbl, 'XTickLabelRotation', 20, ...
	'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10);
colormap(clr_map);
ylabel('Baseline Improvement', 'FontWeight', 'Bold');
legend(Net_lst, 'FontWeight', 'Bold', 'Location', 'NorthWest');

%% Saving
% output_name = sprintf('./Plots/S10_NOPs_Performance.pdf');
% set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [30 10], 'PaperPosition', [0 0 30 10]);
% print('-dpdf', '-r300', output_name);

function Tar_ind = getInd(Src_lst, Tar_lst)
for li=1:numel(Tar_lst)
	Tar_ind(li,1) = find(ismember(Src_lst, Tar_lst{li}));
end
end
