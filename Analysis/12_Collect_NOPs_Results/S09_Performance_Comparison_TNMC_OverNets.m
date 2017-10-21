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
net_cnst = 'P10000';
% net_cnst = 'G00500';
res_info = dir([res_path 'MRK_CVT01_TNMC_*' net_cnst '_MSN-500.mat']);
res_lst = {
    'MRK_CVT01_TNMC_None-G11748_MSN-500.mat'
    'MRK_CVT01_TNMC_HPRD-P10000_MSN-500.mat'
    'MRK_CVT01_TNMC_I2D-P10000_MSN-500.mat'
    'MRK_CVT01_TNMC_MSigDB-P10000_MSN-500.mat'
    'MRK_CVT01_TNMC_STRING-P10000_MSN-500.mat'
    'MRK_CVT01_TNMC_KEGG-P10000_MSN-500.mat'
    'MRK_CVT01_TNMC_AbsCorr-P10000_MSN-500.mat'
    'MRK_CVT01_TNMC_Random-P10000_MSN-500.mat'
    'MRK_CVT01_TNMC_AvgSynACr-P10000_MSN-500.mat'
    }';
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = winter(n_res);
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
    
    dist_h = distributionPlot(auc_rep', 'color', clr_map(si,:), 'xValues', si, 'showMM', 5);
    set(dist_h{2}, 'Color', 'k');
    set(dist_h{2}(1), 'Marker', '.', 'LineWidth', 2, 'MarkerSize', 10);
    %bar(si, mean(auc_rep), 'FaceColor', clr_map(si,:));
    %errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 2, 0.2, [0 0 0]);
    X_lbl{si,1} = sprintf('%s-n%s', res_info{4}, res_info{5}(5:7));
end
xlim([0 n_res+1]);
ylim([0.52 0.67]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', X_lbl, 'XTickLabelRotation', 20, ...
	'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 12, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S09_PerformanceComparison_NMC_OverNet.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperUnits', 'Inches', 'PaperPositionMode','auto', 'PaperSize', [15 5], 'PaperPosition', [0 0 15 5]);
print('-dpdf', '-r300', output_name);

