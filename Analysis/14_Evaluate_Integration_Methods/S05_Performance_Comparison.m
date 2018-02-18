clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

vio_lst = {
%     {'CMB_Rnd_STRING-Shf_NN20_CVT50.mat' 'CMB_Rnd_STRING_NN20_CVT50.mat'}
    {'CMB_Avg_STRING-Shf_NN20_CVT50.mat' 'CMB_Avg_STRING_NN20_CVT50.mat'}
    {'CMB_Std_STRING-Shf_NN20_CVT50.mat' 'CMB_Std_STRING_NN20_CVT50.mat'}
    {'CMB_PCA1_STRING-Shf_NN20_CVT50.mat' 'CMB_PCA1_STRING_NN20_CVT50.mat'}
    {'CMB_DA2_STRING-Shf_NN20_CVT50.mat' 'CMB_DA2_STRING_NN20_CVT50.mat'}
    {'CMB_Reg_STRING-Shf_NN20_CVT50.mat' 'CMB_Reg_STRING_NN20_CVT50.mat'}
%     {'CMB_Reg_STRING-Shf_NN01_CVT51.mat' 'CMB_Reg_STRING_NN01_CVT51.mat'}
%     {'CMB_Reg_STRING-Shf_NN02_CVT51.mat' 'CMB_Reg_STRING_NN02_CVT51.mat'}
%     {'CMB_Reg_STRING-Shf_NN05_CVT51.mat' 'CMB_Reg_STRING_NN05_CVT51.mat'}
%     {'CMB_Reg_STRING-Shf_NN07_CVT51.mat' 'CMB_Reg_STRING_NN07_CVT51.mat'}
%     {'CMB_Reg_STRING-Shf_NN10_CVT51.mat' 'CMB_Reg_STRING_NN10_CVT51.mat'}
%     {'CMB_Reg_STRING-Shf_NN20_CVT51.mat' 'CMB_Reg_STRING_NN20_CVT51.mat'}
%     {'CMB_Reg_STRING-Shf_NN50_CVT51.mat' 'CMB_Reg_STRING_NN50_CVT51.mat'}
};
n_vio = numel(vio_lst);
clr_map = lines(n_vio);
neg_clr = [0.7 0.7 0.7];
y_lim = [0.85 1.15];

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
plot([0 n_vio*3], [1 1], 'r', 'LineWidth', 1.5, 'Color', [246, 82, 82]/255);
hold on
X_lbl = {};
for vi=1:n_vio
    data_l = LoadCmbRes(vio_lst{vi}{1});
    data_r = LoadCmbRes(vio_lst{vi}{2});
    vio_pos(vi) = vi*2;
    
    left_h = distributionPlot(data_l, 'color', neg_clr, 'xValues', vio_pos(vi)-0.4, 'showMM', 0, 'distWidth', 0.7, 'histOri', 'left');
    set(left_h{1}, 'FaceAlpha', 0.8);
    box_h = boxplot(data_l, 'Positions', vio_pos(vi)-0.2, 'Color', neg_clr*0.7, 'Symbol', '', 'Widths', 0.3);
    set(box_h, 'LineWidth', 2);
%     text(vio_pos(vi)-0.2, 1.11, sprintf('n=%d', numel(data_l)), 'HorizontalAlignment', 'Right');
    
    right_h = distributionPlot(data_r, 'color', clr_map(vi,:), 'xValues', vio_pos(vi)+0.4, 'showMM', 0, 'distWidth', 0.7, 'histOri', 'right');
    set(right_h{1}, 'FaceAlpha', 0.8);
    box_h = boxplot(data_r, 'Positions', vio_pos(vi)+0.2, 'Color', clr_map(vi,:)*0.7, 'Symbol', '', 'Widths', 0.3);
    set(box_h, 'LineWidth', 2);
%     text(vio_pos(vi)+0.2, 1.11, sprintf('n=%d', numel(data_l)), 'HorizontalAlignment', 'Left');
    
    [~, pval] = ttest2(data_l, data_r, 'Tail', 'Left');
    
    res_info = regexp(vio_lst{vi}{1}, '_', 'split');
    net_name = sprintf('%s-%s\n%2.0e', res_info{2}, res_info{4}, pval);
    text(vio_pos(vi), y_lim(1), net_name, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', ...
        'FontSize', 12, 'FontWeight', 'Bold');
end
% legend([left_h{1} right_h{1}], {'Random', 'STRING'});
xlim([1 vio_pos(vi)+1]);
ylim(y_lim);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.2f', y), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', vio_pos, 'XTickLabel', [], 'XTickLabelRotation', 0, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.2);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S01_PerformanceComparison.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [16 4], 'PaperPosition', [0 0 16 4]);
print('-dpdf', '-r300', output_name);


function res = LoadCmbRes(cmb_name)
cmb_path = './Combined_AUC/';
cmb_name = [cmb_path cmb_name];
fprintf('Reading [%s]\n', cmb_name);
cmb_data = load(cmb_name);
if any(isnan(cmb_data.Ind_TeAUC(:))) || any(isnan(cmb_data.Cmb_TeAUC(:))), error(); end

mean_ind = median(median(cmb_data.Ind_TeAUC, 3),2);
mean_grp = median(median(cmb_data.Cmb_TeAUC, 3),2);
res = mean_grp./mean_ind;
end