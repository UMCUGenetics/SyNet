clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    'MRK_CVT01_SVM-RBF_None-G11748_MSN-500.mat'
    'MRK_CVT01_SVM-Lin_None-G11748_MSN-500.mat'
    'MRK_CVT01_KNN0_None-G11748_MSN-500.mat'
    'MRK_CVT01_TRnFrst_None-G11748_MSN-500.mat'
    'MRK_CVT01_LDA_None-G11748_MSN-500.mat'
    'MRK_CVT01_TNB_None-G11748_MSN-500.mat'
%     'MRK_CVT01_TReg_None-G11748_MSN-100.mat'
    'MRK_CVT01_TNMC_None-G11748_MSN-050.mat'
%     'MRK_CVT01_TNMC_None-G11748_MSN-020.mat'
    'MRK_CVT01_TReg_None-G11748_MSN-050.mat'
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
};
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(n_res);
clr_map(n_res-2,:) = [0.95 0.95 0];
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
    
    box_h = boxplot(auc_rep', 'Positions', si, 'Width', 0.5, 'Colors', clr_map(si,:), 'Notch', 'off');
    set(box_h, 'LineWidth', 2);
    XData = get(box_h, 'XData');
    YData = get(box_h, 'YData');
    patch(XData{6}([1 2 2 1]), YData{5}([1 1 2 2]), 'b', 'FaceColor', clr_map(si,:), 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    median_lines = findobj(box_h, 'type', 'line', 'Tag', 'Median');
    set(median_lines, 'Color', clr_map(si,:)*0.8);
    uistack(box_h,'top');
    
    X_lbl{si,1} = sprintf('%s-%s', res_info{3}, res_info{5}(5:7));
end
xlim([0 n_res+1]);
ylim([0.52 0.67]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', X_lbl, 'XTickLabelRotation', 10, ...
	'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S02_01_Vanilla_PerformanceComparison.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [10 3], 'PaperPosition', [0 0 10 3]);
print('-dpdf', '-r300', output_name);


