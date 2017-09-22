clc;
clear;
close all

%% Initialization
% addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/Hatchfill/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    'MRK_CVT01_Lasso_Random-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_Random-P10000_MSN-500.mat'
    
    'MRK_CVT01_Lasso_I2D-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_I2D-P10000_MSN-500.mat'
    
    'MRK_CVT01_Lasso_I2D-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_I2D-G11748_MSN-500.mat'
    
    'MRK_CVT01_Lasso_MSigDB-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_MSigDB-P10000_MSN-500.mat'
    
    'MRK_CVT01_Lasso_HPRD-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_HPRD-P10000_MSN-500.mat'
    
    'MRK_CVT01_Lasso_HPRD-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_HPRD-G11748_MSN-500.mat'
    
    'MRK_CVT01_Lasso_KEGG-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_KEGG-P10000_MSN-500.mat'
    
    'MRK_CVT01_Lasso_KEGG-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_KEGG-G11748_MSN-500.mat'
    
    'MRK_CVT01_Lasso_STRING-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_STRING-P10000_MSN-500.mat'
    
    'MRK_CVT01_Lasso_STRING-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_STRING-G11748_MSN-500.mat'
    
    'MRK_CVT01_Lasso_ACr-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_ACr-P10000_MSN-500.mat'
    
    'MRK_CVT01_Lasso_AvgSynACr-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_AvgSynACr-P10000_MSN-500.mat'
    
%     'MRK_CVT01_GLasso_AvgSynACr-P10000_MSN-500.mat'
%     'MRK_CVT01_GLasso10_AvgSynACr-P10000_MSN-500.mat'
    
%     'MRK_CVT01_GLasso2_AvgSynACr-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_AvgSynACr-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso7_AvgSynACr-P10000_MSN-500.mat'
%     'MRK_CVT01_GLasso7_AvgSynACr-G05000_MSN-500.mat'
%     'MRK_CVT01_GLasso10_AvgSynACr-P10000_MSN-500.mat'
%     'MRK_CVT01_GLasso20_AvgSynACr-P10000_MSN-500.mat'
    
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
    };
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(ceil(n_res/2));
X_lbl = {};
AUC_cmb = [];
xlim([0 n_res+1]);
ylim([0.60 0.68]);
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
    %bar_h = bar(si, mean(auc_rep), 'FaceColor', clr_map(ceil(si/2),:));
    bar_h = patch([si-0.25 si+0.25 si+0.25 si-0.25], [0 0 mean(auc_rep) mean(auc_rep)], 'b', 'FaceColor', clr_map(ceil(si/2),:)); 
    if mod(si,2)==0
        patch_h = findobj(bar_h, 'Type', 'patch');
        hatch_h = hatchfill(patch_h, 'single', -45, 5, clr_map(ceil(si/2),:));
        hatch_h.Color = [0 0 0];
    end
    errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 1.5, 0.1, [0 0 0]);
    net_info = regexp(res_info{4}, '-', 'Split');
    X_lbl{si,1} = sprintf('%s', net_info{1});
end
text(n_res/2+0.5, 0.68, {'Sparse Group Lasso' '(10000 pairs)'}, 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', X_lbl, 'XTickLabelRotation', 10, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
% output_name = sprintf('./Plots/S07_PerformanceComparison_NetBased_KeepPairsConstant.pdf');
% set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [13 3], 'PaperPosition', [0 0 13 3]);
% print('-dpdf', '-r300', output_name);

