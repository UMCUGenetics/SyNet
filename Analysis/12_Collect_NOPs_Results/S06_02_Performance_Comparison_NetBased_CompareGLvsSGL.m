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
    'MRK_CVT01_GLasso_Random-G00500_MSN-500.mat'
    'MRK_CVT01_CFGLasso_Random-G00500_MSN-500.mat'
    
    'MRK_CVT01_GLasso_I2D-G00500_MSN-500.mat'
    'MRK_CVT01_CFGLasso_I2D-G00500_MSN-500.mat'
    
    'MRK_CVT01_GLasso_MSigDB-G00500_MSN-500.mat'
    'MRK_CVT01_CFGLasso_MSigDB-G00500_MSN-500.mat'
    
    'MRK_CVT01_GLasso_HPRD-G00500_MSN-500.mat'
    'MRK_CVT01_CFGLasso_HPRD-G00500_MSN-500.mat'
    
    'MRK_CVT01_GLasso_KEGG-G00500_MSN-500.mat'
    'MRK_CVT01_CFGLasso_KEGG-G00500_MSN-500.mat'
    
    'MRK_CVT01_GLasso_STRING-G00500_MSN-500.mat'
    'MRK_CVT01_CFGLasso_STRING-G00500_MSN-500.mat'
    
    'MRK_CVT01_GLasso_ACr-G00500_MSN-500.mat'
    'MRK_CVT01_CFGLasso_ACr-G00500_MSN-500.mat'
    
    'MRK_CVT01_GLasso_AvgSynACr-G00500_MSN-500.mat'
    'MRK_CVT01_CFGLasso_AvgSynACr-G00500_MSN-500.mat'
    
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
y_lim = [0.60 0.68];
ylim(y_lim);
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
    
    net_info = regexp(res_info{4}, '-', 'Split');
    if isequal(res_info{3}, 'CFGLasso') %mod(si,2)==0
        patch_h = findobj(bar_h, 'Type', 'patch');
        hatch_h = hatchfill(patch_h, 'single', -45, 5, clr_map(ceil(si/2),:));
        X_lbl{si,1} = sprintf('%s\nSGL', net_info{1});
    else
        X_lbl{si,1} = sprintf('%s\nGL', net_info{1});
    end
    errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 1.5, 0.1, [0 0 0]);
end
X_lbl(end-2:end) = {sprintf('SyNet\nGL') sprintf('SyNet\nSGL') 'All genes'};
for si=1:n_res
    text(si, y_lim(1), X_lbl{si}, 'FontSize', 9, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
end
% text(n_res/2+0.5, 0.68, 'Sparse Group Lasso', 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
% legend([patch_h hatch_h], {'Group Lasso' 'Sparse Group Lasso'});

y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', [], ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S06_02_PerformanceComparison_NetBased_GL_vs_SGL.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [13 3], 'PaperPosition', [0 0 13 3]);
print('-dpdf', '-r300', output_name);
