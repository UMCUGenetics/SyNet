clc;
clear;
close all

%% Initialization
% addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/Hatchfill/');
addpath('../_Utilities/');

%% Select methods and networks
res_lst = {
    {'MRK_CVT01_iPark_AvgSynACr-P50000_MSN-500.mat' 'MRK_CVT01_SFiPark_AvgSynACr-P50000_MSN-500'}    
    {'MRK_CVT01_iChuang_AvgSynACr-P50000_MSN-500.mat' 'MRK_CVT01_SFiChuang_AvgSynACr-P50000_MSN-500'}
    {'MRK_CVT01_iTaylor_AvgSynACr-P50000_MSN-500.mat' 'MRK_CVT01_SFiTaylor_AvgSynACr-P50000_MSN-500'}
    {'MRK_CVT01_NetGL_AvgSynACr-P50000_MSN-500.mat' 'MRK_CVT01_NetSFGL_AvgSynACr-P50000_MSN-500.mat'}
    {'MRK_CVT01_LExAG_None-G11748_MSN-500.mat' 'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'}
    };
n_res = numel(res_lst);
y_lim = [0.60 0.67];

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
X_lbl = {};
clr_map = [
    0.3490    0.6980    0.8706
    0.8863    0.6039    0.4157
    0.8706    0.3490    0.6549
    1.0000    0.1000    0.1000
    0.4000    0.6000    0.4000
    ];

xlim([0.5 n_res+0.5]);
ylim(y_lim);
for si=1:n_res
    auc1_study = LoadData(res_lst{si}{1});
    auc2_study = LoadData(res_lst{si}{2});
    res_info = regexp(res_lst{si}{2}, '_', 'split');
    net_info = regexp(res_info{4}, '-', 'Split');
    X_lbl{si,1} = sprintf('%s', net_info{1});
    met_clr = clr_map(si,:);
    
    bar_pos = si-0.2;
    bar_h = patch(bar_pos+[-1 1 1 -1]*0.18, [0 0 mean(auc1_study) mean(auc1_study)], 'b', 'FaceColor', met_clr);
    errorbarEx(bar_pos, mean(auc1_study), std(auc1_study), std(auc1_study), 1.5, 0.1, [0 0 0]);
    %text(bar_pos, mean(auc1_study), sprintf('%0.2f', std(auc1_study)*100), 'FontSize', 12, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');
    
    if si~=n_res
        bar_pos = si+0.2;
        bar_h = patch(bar_pos+[-1 1 1 -1]*0.18, [0 0 mean(auc2_study) mean(auc2_study)], 'b', 'FaceColor', met_clr);
        patch_h = findobj(bar_h, 'Type', 'patch');
        hatch_h = hatchfill(patch_h, 'single', -45, 10, met_clr*0.6+ ones(1,3)*0.2);
        set(hatch_h, 'Color', [0.9 0.9 0.9], 'LineWidth', 1.5);
        errorbarEx(bar_pos, mean(auc2_study), std(auc2_study), std(auc2_study), 1.5, 0.1, [0 0 0]);
        %text(bar_pos, mean(auc2_study), sprintf('%0.2f', std(auc2_study)*100), 'FontSize', 12, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Left', 'VerticalAlignment', 'Bottom');
    end
    text(si, y_lim(1)-0.002, res_info{3}, 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top', 'Rotation', 0);
end
legend([bar_h, patch_h], {'Network based', 'Shuffled nodes'});

y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', [], 'XTickLabel', [], 'XTickLabelRotation', 10, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 7, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

% return
%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_06_NetVsShuffle.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [7 3], 'PaperPosition', [0 0 7 3]);
print('-dpdf', '-r300', output_name);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function auc_study = LoadData(res_name)
res_path = './Collected_Results/';

res_name = [res_path res_name];
fprintf('Reading [%s]\n', res_name);
res_data = load(res_name);
if any(isnan(res_data.AUC_mat(:))), error(); end
auc_study = mean(res_data.AUC_mat,1);
end
