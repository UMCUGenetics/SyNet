clc;
clear;
close all

%% Initialization
% addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    'MRK_CVT01_TLEx_Random-G00700_MSN-700.mat'
    'MRK_CVT01_TLEx_KEGG-G11748_MSN-700.mat'
    'MRK_CVT01_TLEx_MSigDB-G11748_MSN-700.mat'
    'MRK_CVT01_TLEx_STRING-G11748_MSN-700.mat'
    'MRK_CVT01_TLEx_ACr-G11748_MSN-700.mat'
    'MRK_CVT01_TLEx_HPRD-G11748_MSN-700.mat'
    'MRK_CVT01_TLEx_I2D-G11748_MSN-700.mat'
    'MRK_CVT01_TLEx_None-G11748_MSN-700.mat'
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
    'MRK_CVT01_TLEx_AvgSynACr-P10000_MSN-700.mat'
    };
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = min(hsv(n_res)*1.0, ones(n_res,3));
X_lbl = {};
AUC_cmb = [];
for si=1:n_res
    res_name = [res_path res_lst{si}];
    res_info = regexp(res_lst{si}, '_', 'split');
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    auc_rep = mean(res_data.AUC_mat);
    AUC_cmb(:,si) = auc_rep;
    
    %dist_h = distributionPlot(auc_rep', 'color', clr_map(si,:), 'xValues', si, 'showMM', 5);
    %set(dist_h{2}, 'Color', 'k');
    %set(dist_h{2}(1), 'Marker', '.', 'LineWidth', 2, 'MarkerSize', 10);
    bar(si, mean(auc_rep), 'FaceColor', clr_map(si,:), 'BarWidth', 0.7, 'EdgeColor', clr_map(si,:)*0.6);
    errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 2, 0.1, clr_map(si,:)*0.5);
    net_info = regexp(res_info{4}, '-', 'Split');
    X_lbl{si,1} = sprintf('%s', net_info{1});
end
X_lbl(end-2:end) = {'ttest-#700' 'All genes' 'SyNet'};
% text(n_res/2+0.5, 0.66, 'Lasso with all genes', 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
xlim([0 n_res+1]);
ylim([0.6 0.66]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', X_lbl, 'XTickLabelRotation', 10, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S05_PerformanceComparison_GenesInNet_G700.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [13 3], 'PaperPosition', [0 0 13 3]);
print('-dpdf', '-r300', output_name);

