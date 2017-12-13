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
    'MRK_CVT01_NetLasso_BioPlex-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_BioPlex-P25000_MSN-500.mat'
    
    'MRK_CVT01_NetLasso_BioGRID-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_BioGRID-P25000_MSN-500.mat'
    
    'MRK_CVT01_NetLasso_IntAct-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_IntAct-P25000_MSN-500.mat'
    
    'MRK_CVT01_NetLasso_STRING-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_STRING-P25000_MSN-500.mat'
        
    'MRK_CVT01_NetLasso_HBGland-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_HBGland-P25000_MSN-500.mat'

    'MRK_CVT01_NetLasso_HBLympNode-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_HBLympNode-P25000_MSN-500.mat'
    
    'MRK_CVT01_NetLasso_HBEpith-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_HBEpith-P25000_MSN-500.mat'
        
    'MRK_CVT01_NetLasso_ACr-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_ACr-P25000_MSN-500.mat'
        
    'MRK_CVT01_NetLasso_Syn-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_Syn-P25000_MSN-500.mat'
    
    'MRK_CVT01_NetLasso_AvgSyn-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_AvgSyn-P25000_MSN-500.mat'
    
    'MRK_CVT01_NetLasso_AvgSynACr-P25000_MSN-500.mat'
    'MRK_CVT01_NetGL_AvgSynACr-P25000_MSN-500.mat'

%     'MRK_CVT01_Lasso_AvgSynACr-P10000_MSN-500.mat'
%     'MRK_CVT01_GLasso_AvgSynACr-P10000_MSN-500.mat'
    
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
    };
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(ceil(n_res/2));
X_lbl = {};
xlim([0 n_res+1]);
ylim([0.61 0.67]);
for si=1:n_res
    res_name = [res_path res_lst{si}];
    res_info = regexp(res_lst{si}, '_', 'split');
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    auc_rep = mean(res_data.AUC_mat,1);
    
    bar_h = patch([si-0.25 si+0.25 si+0.25 si-0.25], [0 0 mean(auc_rep) mean(auc_rep)], 'b', 'FaceColor', clr_map(ceil(si/2),:)); 
    if isequal(res_info{3}, 'NetGL')
        patch_h = findobj(bar_h, 'Type', 'patch');
        hatchfill(patch_h, 'single', -45, 5, clr_map(ceil(si/2),:));
    end
    errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 1.5, 0.1, [0 0 0]);
    net_info = regexp(res_info{4}, '-', 'Split');
    X_lbl{si,1} = sprintf('%s', net_info{1});
end
X_lbl(end-2:end) = {'SyNet' 'SyNet' 'All genes'};
% text(n_res/2+0.5, 0.68, 'Sparse Group Lasso', 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');

y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', X_lbl, 'XTickLabelRotation', 10, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

return
%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_NetBased.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [13 3], 'PaperPosition', [0 0 13 3]);
print('-dpdf', '-r300', output_name);

