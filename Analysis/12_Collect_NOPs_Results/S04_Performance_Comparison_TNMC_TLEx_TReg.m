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
    'MRK_CVT01_TNMC_None-G11748_MSN-020.mat'
    'MRK_CVT01_TReg_None-G11748_MSN-020.mat'
    'MRK_CVT01_TLEx_None-G11748_MSN-020.mat'
    'MRK_CVT01_TNMC_None-G11748_MSN-050.mat'
    'MRK_CVT01_TReg_None-G11748_MSN-050.mat'
    'MRK_CVT01_TLEx_None-G11748_MSN-050.mat'
    'MRK_CVT01_TNMC_None-G11748_MSN-100.mat'
    'MRK_CVT01_TReg_None-G11748_MSN-100.mat'
    'MRK_CVT01_TLEx_None-G11748_MSN-100.mat'
    'MRK_CVT01_TNMC_None-G11748_MSN-500.mat'
    'MRK_CVT01_TReg_None-G11748_MSN-500.mat'
    'MRK_CVT01_TLEx_None-G11748_MSN-500.mat'
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
};
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = [
    0.3000    0.7470    1.0000
    0.3333    1.0000    0.6667
    0.6706    0.4118    0.9686
    0.9706    0.6118    0.4686
    0.9506    0.3018    0.3018
    ];
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
    
    dist_h = distributionPlot(auc_rep', 'color', clr_map(ceil(si/3),:), 'xValues', si, 'showMM', 5);
    set(dist_h{2}, 'Color', 'k');
    set(dist_h{2}(1), 'Marker', '.', 'LineWidth', 2, 'MarkerSize', 10);
    %bar(si, mean(auc_rep), 'FaceColor', clr_map(si,:));
    %errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 2, 0.2, [0 0 0]);
    X_lbl{si,1} = sprintf('%s-%s', res_info{3}, res_info{5}(5:7));
    if mod(si,3)==2
        text(si, 0.66, sprintf('#Gene=%s', res_info{5}(5:7)), 'FontSize', 12, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center');
    end
end
xlim([0 n_res+1]);
ylim([0.52 0.67]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', X_lbl, 'XTickLabelRotation', 20, ...
	'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S04_PerformanceComparison_NMC-TLEx-TReg.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [10 3], 'PaperPosition', [0 0 10 3]);
print('-dpdf', '-r300', output_name);

