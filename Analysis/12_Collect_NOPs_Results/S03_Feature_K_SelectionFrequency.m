clc;
clear;
close all

%% Initialization
% addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Distribution Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Load data
res_path = './Collected_Results/';
nmc_data = load([res_path 'MRK_CVT01_TNMCAd_None-G11748_MSN-500.mat']);
las_data = load([res_path 'MRK_CVT01_LExAG_None-G11748_MSN-500.mat']);
knn_data = load([res_path 'MRK_CVT01_KNN0_None-G11748_MSN-500.mat']);

%% Plotting NMC and Lasso
close all
figure('Position', [100 100 1500 400]);

subplot(1,2,1);
hold on
Opt_K = sort(nmc_data.Opt_K);
dist_h = distributionPlot(Opt_K(1:100), 'color', [192, 134, 253]/255, 'xValues', 1, 'showMM', 5);
set(dist_h{2}, 'Color', 'k');
delete(dist_h{2}(1));
plot(1, median(Opt_K), 'Marker', 'd', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

n_nonzero = sort(sum(las_data.Marker_lst>0, 2));
dist_h = distributionPlot(n_nonzero(1:130), 'color', [69, 252, 90]/255, 'xValues', 2, 'showMM', 5);
set(dist_h{2}, 'Color', 'k');
delete(dist_h{2}(1));
plot(2, median(n_nonzero), 'Marker', 'd', 'LineWidth', 2, 'MarkerSize', 7, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');

xlim([0.5 2.5]);
ylim([0 120]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%d', y), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:2, 'XTickLabel', {'NMC' 'Lasso'}, 'XTickLabelRotation', 0, ...
	'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 12, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.2);
ylabel('# Features selected', 'FontWeight', 'Bold');

subplot(1,2,2);
hold on
[lst, freq] = getTop(knn_data.Opt_K);
n_elm = numel(lst);
for i=1:n_elm
    bar(i, freq(i), 'FaceColor', [252, 142, 69]/255, 'BarWidth', 0.8);
end
xlim([0 n_elm+1]);
set(gca, 'XTick', 1:n_elm, 'XTickLabel', lst, 'XTickLabelRotation', 0, ...
	'FontWeight', 'Bold', 'FontSize', 12, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.2);
ylabel('Frequency of selected K', 'FontWeight', 'Bold');
title('KNN', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S03_Feature_K_SelectionFrequency.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [15 4], 'PaperPosition', [0 0 15 4]);
print('-dpdf', '-r300', output_name);
