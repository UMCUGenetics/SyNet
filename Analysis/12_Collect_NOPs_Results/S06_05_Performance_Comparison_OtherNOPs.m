clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/Hatchfill/');
addpath('../_Utilities/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Select methods and networks
grp_lst = {
    {
    'MRK_CVT01_iPark_HumanInt-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_BioPlex-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_BioGRID-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_IntAct-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_STRING-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_HBOvary-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_HBBrain-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_HBKidney-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_HBGland-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_HBLympNode-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_ACr-P50000_MSN-500.mat'
    'MRK_CVT01_iPark_AvgSynACr-P10000_MSN-500.mat'}
    {
    'MRK_CVT01_iChuang_HumanInt-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_BioPlex-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_BioGRID-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_IntAct-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_STRING-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_HBOvary-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_HBBrain-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_HBKidney-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_HBGland-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_HBLympNode-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_ACr-P50000_MSN-500.mat'
    'MRK_CVT01_iChuang_AvgSynACr-P10000_MSN-500.mat'}
    {
    'MRK_CVT01_iTaylor_HumanInt-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_BioPlex-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_BioGRID-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_IntAct-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_STRING-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_HBOvary-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_HBBrain-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_HBKidney-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_HBGland-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_HBLympNode-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_ACr-P50000_MSN-500.mat'
    'MRK_CVT01_iTaylor_AvgSynACr-P10000_MSN-500.mat'}
%     {
%     'MRK_CVT01_NetLasso_HumanInt-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_BioPlex-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_BioGRID-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_IntAct-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_STRING-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_HBOvary-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_HBBrain-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_HBKidney-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_HBGland-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_HBLympNode-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_ACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetLasso_AvgSynACr-P10000_MSN-500.mat'}
    {
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'}
%     {
%     'MRK_CVT01_NetLasso_AvgSynACr-P50000_MSN-500.mat'
%     'MRK_CVT01_NetGL_AvgSynACr-P50000_MSN-500.mat'}
    };
n_grp = numel(grp_lst);
y_lim = [0.60 0.67];
offset = 0.5;

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
X_lbl = {};
ylim(y_lim);
bar_pos = 1;
for gi=1:n_grp
    res_lst = grp_lst{gi};
    n_res = numel(res_lst);
    for ri=1:n_res
        auc_study = LoadData(res_lst{ri});
        res_info = regexp(res_lst{ri}, '_', 'split');
        net_info = regexp(res_info{4}, '-', 'Split');
        met_lbl = regexprep(res_info{3}, '^i', '');
        
        [net_clr, net_lbl{ri}] = getColor(net_info{1});
        bar_h(ri) = patch(bar_pos+[-1 1 1 -1]*0.2, [0 0 mean(auc_study) mean(auc_study)], 'b', ...
            'FaceColor', net_clr, 'EdgeColor', net_clr*0.5);
        % [~, bar_h(ri)] = BoxPlotEx(auc_study, 'Color', net_clr, 'Positions', bar_pos, 'Widths', 0.3, 'Symbol', '');
            
        errorbarEx(bar_pos, mean(auc_study), std(auc_study), std(auc_study), 1.5, 0.09, net_clr*0.6);
        text(bar_pos, y_lim(1)-0.002, net_lbl{ri}, 'FontSize', 8, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Right', ...
            'VerticalAlignment', 'Middle', 'Rotation', 30);
        bar_pos = bar_pos + offset;
    end
    text(bar_pos - 3, 0.65, met_lbl, 'FontSize', 16, 'FontWeight', 'Bold', 'HorizontalAlignment', ...
            'Center', 'VerticalAlignment', 'Bottom', 'Rotation', 0);
    bar_pos = bar_pos + offset + 0.2;
end
legend(bar_h, net_lbl, 'Location', 'NorthWest');
xlim([0.5 bar_pos]);

y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', [], 'XTickLabel', [], 'XTickLabelRotation', 10, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 8, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_05_OtherNOPs.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [15 4], 'PaperPosition', [0 0 15 4]);
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
