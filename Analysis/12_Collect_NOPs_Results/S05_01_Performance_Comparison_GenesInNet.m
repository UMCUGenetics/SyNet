clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../_Utilities/');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    'MRK_CVT01_NetLasso_HumanInt-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_BioPlex-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_BioGRID-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_IntAct-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_STRING-P50000_MSN-500.mat'
    ''
    'MRK_CVT01_NetLasso_HBOvary-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBBrain-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBKidney-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBGland-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBLympNode-P50000_MSN-500.mat'
    ''
    'MRK_CVT01_NetLasso_ACr-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_AvgSynACr-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_AvgSynACrNShuff-P50000_MSN-500.mat'
    ''
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
    };
n_res = numel(res_lst);
y_lim = [0.6 0.66];

%% Plotting performance
close all
figure('Position', [100 100 700 400]);
hold on
x_tick = 1:n_res;
for si=1:n_res
    if isempty(res_lst{si})
        x_tick(si) = nan;
        continue;
    end
    
    res_name = [res_path res_lst{si}];
    res_info = regexp(res_lst{si}, '_', 'split');
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    auc_rep = mean(res_data.AUC_mat, 1);
    net_info = regexp(res_info{4}, '-', 'Split');
    
    [met_clr, Method_lbl{si,1}] = getColor(net_info{1});
    bar(si, mean(auc_rep), 'FaceColor', met_clr, 'BarWidth', 0.7, 'EdgeColor', met_clr*0.6);
    err_h = errorbarEx(si, mean(auc_rep), std(auc_rep), 0, 2, 0.1, met_clr*0.8);
    delete(err_h(3));
    
    if si==n_res
        x_lbl = sprintf('%d\n%0.0f', res_data.Gene_Map.Count, round(median(sum(res_data.Marker_lst~=0,2))));
    else
        x_lbl = sprintf('%d\n%0.0f', mode(res_data.BestNetwork), round(median(sum(res_data.Marker_lst~=0,2))));
    end
    text(si, y_lim(1), x_lbl, 'FontSize', 7, 'FontWeight', 'Bold', 'Color', [0 0 0], ...
        'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Center');
    text(si, y_lim(1), Method_lbl{si}, 'FontSize', 7, 'FontWeight', 'Bold', 'Rotation', 25, ...
        'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Right');
end
x_tick(isnan(x_tick)) = [];

xlim([0 n_res+1]);
ylim(y_lim);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', x_tick, 'XTickLabel', [], ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 7, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

% return
%% Saving
output_name = sprintf('./Plots/S05_01_PerformanceComparison_GenesInNet.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [7 3], 'PaperPosition', [0 0 7 3]);
print('-dpdf', '-r300', output_name);

