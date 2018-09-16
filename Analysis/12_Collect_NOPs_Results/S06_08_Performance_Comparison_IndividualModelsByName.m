clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/Hatchfill/');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
%     'MRK_SyN10-SyN10_CVT01_LExAG_None-G11748_MSN-500.mat'
    'MRK_SyN25-SyN25_CVT01_LExAG_None-G11748_MSN-500.mat'
    'MRK_SyN50-SyN50_CVT01_LExAG_None-G11748_MSN-500.mat'
%     'MRK_SyN100-SyN100_CVT01_LExAG_None-G11748_MSN-500.mat'
    'MRK_SyNet-SyNet_CVT01_LExAG_None-G11748_MSN-500.mat'
%     'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
%     ''
%     'MRK_CVT01_HubLasso5_SyNet-AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_HubGL5_SyNet-AvgSynACr-P50000_MSN-500.mat'
%     ''
%     'MRK_SyN10-SyN10_CVT01_NetLasso_STRING-P50000_MSN-500.mat'
%     'MRK_SyN25-SyN25_CVT01_NetLasso_STRING-P50000_MSN-500.mat'
%     'MRK_SyN50-SyN50_CVT01_NetLasso_STRING-P50000_MSN-500.mat'
%     'MRK_SyN100-SyN100_CVT01_NetLasso_STRING-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_STRING-P50000_MSN-500.mat'
    ''
%     'MRK_SyN10-SyN10_CVT01_NetGL_STRING-P50000_MSN-500.mat'
    'MRK_SyN25-SyN25_CVT01_NetGL_STRING-P50000_MSN-500.mat'
    'MRK_SyN50-SyN50_CVT01_NetGL_STRING-P50000_MSN-500.mat'
%     'MRK_SyN100-SyN100_CVT01_NetGL_STRING-P50000_MSN-500.mat'
    'MRK_SyNet-SyNet_CVT01_NetGL_STRING-P50000_MSN-500.mat'
%     'MRK_CVT01_NetGL_STRING-P50000_MSN-500.mat'
%     ''
%     'MRK_CVT01_NetLasso_SyN10-ACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_SyN25-ACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_SyN50-ACr-P50000_MSN-500.mat'
    ''
%     'MRK_CVT01_NetGL_SyN10-ACr-P50000_MSN-500.mat'
    'MRK_SyN25-SyN25_CVT01_NetGL_SyN25-ACr-P50000_MSN-500.mat'
    'MRK_SyN50-SyN50_CVT01_NetGL_SyN50-ACr-P50000_MSN-500.mat'
    'MRK_SyNet-SyNet_CVT01_NetGL_ACr-P50000_MSN-500.mat'
%     ''
%     'MRK_CVT01_NetLasso_SyN10-AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_SyN25-AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_SyN50-AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_SyN100-AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_SyNet-AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_AvgSynACr-P50000_MSN-500.mat'
    ''
%     'MRK_CVT01_NetGL_SyN10-AvgSynACr-P50000_MSN-500.mat'
    'MRK_SyN25-SyN25_CVT01_NetGL_SyN25-AvgSynACr-P50000_MSN-500.mat'
    'MRK_SyN50-SyN50_CVT01_NetGL_SyN50-AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetGL_SyN100-AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetGL_SyNet-AvgSynACr-P50000_MSN-500.mat'
    'MRK_SyNet-SyNet_CVT01_NetGL_AvgSynACr-P50000_MSN-500.mat'
    };
n_res = numel(res_lst);
offset = 0.4;

%% Plotting performance
close all
figure('Position', [100 100 800 400]);
hold on
X_lbl = {};
AUC_cmb = [];
y_lim = [0.60 0.67];
xlim([0 n_res+1]);
ylim(y_lim);
for si=1:n_res
    if isempty(res_lst{si})
        continue;
    end
    res_name = [res_path res_lst{si}];
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    auc_rep = mean(res_data.AUC_mat, 1);
    res_info = regexp(res_lst{si}, '_', 'split');
    net_info = regexp(res_info{4}, '-', 'Split');
    [met_clr, met_lbl] = getColor(res_info{5});
    
    if mod(si, 4)==2
        text(si, 0.67, sprintf('%s\n%s', res_info{2}, res_info{5}), 'FontWeight', 'Bold');
    end
    
    bar_h = patch(si+[-1 1 1 -1]*offset, [0 0 mean(auc_rep) mean(auc_rep)], 'b', 'FaceColor', met_clr);
    %     if ~strcmp(res_info{3}, 'Lasso')
    %         patch_h = findobj(bar_h, 'Type', 'patch');
    %         hatch_h = hatchfill(patch_h, 'single', -45, 5, met_clr);
    %         hatch_h.Color = [0 0 0];
    %     end
    errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 1.5, 0.1, [0 0 0]);
    
    text(si, y_lim(1)*0.999, sprintf('%s\n%s', res_info{2}, met_lbl), 'FontSize', 12, 'FontWeight', 'Bold', 'Rotation', 0, ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
end
%text(n_res/2+0.5, 0.68, {'Group Lasso' '(10000 pairs)'}, 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');

y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', [], 'XTickLabel', [], 'XTickLabelRotation', 0, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

return
%% Saving
output_name = sprintf('./Plots/S07_PerformanceComparison_08_ByModelName.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [12 4], 'PaperPosition', [0 0 12 4]);
print('-dpdf', '-r300', output_name);

