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
cmb_path = './Combined_AUC/';
cmb_lst = {
    'CMB_Rnd_Random_NN05_CVT50.mat'
    'CMB_Rnd_STRING_NN05_CVT50.mat'
    'CMB_Avg_Random_NN05_CVT50.mat'
    'CMB_Avg_STRING_NN05_CVT50.mat'
    'CMB_Std_Random_NN05_CVT50.mat'
    'CMB_Std_STRING_NN05_CVT50.mat'
%     'CMB_DA2NoRem_Random_NN05_CVT50.mat'
%     'CMB_DA2NoRem_STRING_NN05_CVT50.mat'
%     'CMB_DPCA_Random_NN05_CVT50.mat'
%     'CMB_DPCA_STRING_NN05_CVT50.mat'
    'CMB_PCA1_Random_NN05_CVT50.mat'
    'CMB_PCA1_STRING_NN05_CVT50.mat'
    'CMB_DA2_Random_NN05_CVT50.mat'
    'CMB_DA2_STRING_NN05_CVT50.mat'
    'CMB_Reg_Random_NN05_CVT50.mat'
    'CMB_Reg_STRING_NN05_CVT50.mat'
    };
n_cmb = numel(cmb_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(n_cmb/2);
X_lbl = {};
for ci=1:n_cmb
    clr_ind = ceil(ci/2);
    cmb_name = [cmb_path cmb_lst{ci}];
    res_info = regexp(cmb_lst{ci}, '_', 'split');
    fprintf('Reading [%s]\n', cmb_name);
    cmb_data = load(cmb_name);
    if any(isnan(cmb_data.Ind_AUC(:))) || any(isnan(cmb_data.Te_AUC(:))), error(); end
    mean_ind = mean(mean(cmb_data.Ind_AUC, 3),2);
    mean_grp = mean(mean(cmb_data.Te_AUC, 3),2);
    IMP_cmb{ci} = mean_grp./mean_ind;
    
    %top_imp = IMP_cmb(:,ci);
    %top_imp = sort(IMP_cmb(:,ci), 'Descend');
    
    dist_h = distributionPlot(IMP_cmb{ci}, 'color', clr_map(clr_ind,:), 'xValues', ci, 'showMM', 5);
    plot(ci, max(IMP_cmb{ci}), 'bx', 'MarkerSize', 12, 'MarkerEdgeColor', clr_map(clr_ind,:)*0.7, 'LineWidth', 2);
    set(dist_h{2}, 'Color', 'k');
    set(dist_h{2}(1), 'Marker', '.', 'LineWidth', 2, 'MarkerSize', 12);
    %bar(si, mean(auc_rep), 'FaceColor', clr_map(si,:));
    %errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 2, 0.2, [0 0 0]);
    X_lbl{ci,1} = sprintf('%s\nn=%d', res_info{3}, numel(IMP_cmb{ci})); % , res_info{3}
    text(ci, 0.9, X_lbl{ci}, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', ...
            'FontSize', 10, 'FontWeight', 'Bold');
    if mod(ci,2)==0
        [~, pval] = ttest2(IMP_cmb{ci-1}, IMP_cmb{ci}, 'Tail', 'Left');
        net_name = sprintf('%s\n%2.0e', res_info{2}, pval);
        text(ci-0.5, 1.11, net_name, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', ...
            'FontSize', 12, 'FontWeight', 'Bold');
    end
%     text(ci-0.1, 1.05, sprintf('n=%d', numel(IMP_cmb{ci})), 'HorizontalAlignment', 'Right', ...
%             'FontSize', 10, 'FontWeight', 'Bold');
end
xlim([0 n_cmb+1]);
ylim([0.90 1.13]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.2f', y), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_cmb, 'XTickLabel', [], 'XTickLabelRotation', 0, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');
plot(xlim, [1 1], 'r', 'LineWidth', 1.5, 'Color', [246, 82, 255]/255);

%% Saving
output_name = sprintf('./Plots/S01_PerformanceComparison.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [16 4], 'PaperPosition', [0 0 16 4]);
print('-dpdf', '-r300', output_name);

