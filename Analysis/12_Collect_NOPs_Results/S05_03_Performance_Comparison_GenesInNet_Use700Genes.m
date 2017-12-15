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
    'MRK_CVT01_NetLasso_BioPlex-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_BioGRID-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_IntAct-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_STRING-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_STRINGnShuff-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBBrain-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBGland-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBLympNode-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_ACr-P25000_MSN-500.mat'
    'MRK_CVT01_NetLasso_AvgSynACr-P25000_MSN-500.mat'
    'MRK_CVT01_TLEx_tTest-G11748_MSN-700.mat'
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
    };
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
for si=1:n_res
    res_name = [res_path res_lst{si}];
    res_info = regexp(res_lst{si}, '_', 'split');
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    auc_rep = mean(res_data.AUC_mat, 1);
    net_info = regexp(res_info{4}, '-', 'Split');
    X_lbl = sprintf('%s', net_info{1});
    
    [met_clr, Method_lbl{si,1}] = getColor(X_lbl);
    bar(si, mean(auc_rep), 'FaceColor', met_clr, 'BarWidth', 0.7, 'EdgeColor', met_clr*0.6);
    errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 2, 0.1, met_clr*0.5);
end

xlim([0 n_res+1]);
ylim([0.6 0.66]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', Method_lbl, 'XTickLabelRotation', 10, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

return
%% Saving
output_name = sprintf('./Plots/S05_PerformanceComparison_GenesInNet_G700.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [13 3], 'PaperPosition', [0 0 13 3]);
print('-dpdf', '-r300', output_name);

