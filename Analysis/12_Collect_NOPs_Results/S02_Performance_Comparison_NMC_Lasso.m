clc;
clear;
close all

%% Initialization
% addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Distribution Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    'MRK_CVT01_TNMC_None-G11748_MSN-020.mat'
    'MRK_CVT01_TNMC_None-G11748_MSN-050.mat'
    'MRK_CVT01_TNMC_None-G11748_MSN-100.mat'
    'MRK_CVT01_TNMC_None-G11748_MSN-500.mat'
    'MRK_CVT01_TNMCAd_None-G11748_MSN-500.mat'
    'MRK_CVT01_KNN1_None-G11748_MSN-500.mat'
    'MRK_CVT01_KNN3_None-G11748_MSN-500.mat'
    'MRK_CVT01_KNN5_None-G11748_MSN-500.mat'
    'MRK_CVT01_KNN7_None-G11748_MSN-500.mat'
    'MRK_CVT01_KNN0_None-G11748_MSN-500.mat'
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
%     'MRK_CVT01_CFGLasso_AvgSynACr-P10000_MSN-500.mat'
%     'MRK_CVT01_CFGLasso_AvgSynACr-G00500_MSN-500.mat'
%     'MRK_CVT01_TLEx_AvgSynACr-P10000_MSN-500.mat'
};
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(n_res);
X_lbl = {};
AUC_cmb = [];
for si=1:n_res
    res_name = [res_path res_lst{si}];
    res_info = regexp(res_lst{si}, '_', 'split');
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    auc_rep = mean(res_data.AUC_mat);
    AUC_cmb(:,si) = auc_rep;
    
    dist_h = distributionPlot(auc_rep', 'color', clr_map(si,:), 'xValues', si, 'showMM', 5);
    set(dist_h{2}, 'Color', 'k');
    set(dist_h{2}(1), 'Marker', '.', 'LineWidth', 2, 'MarkerSize', 10);
    %bar(si, mean(auc_rep), 'FaceColor', clr_map(si,:));
    %errorbarEx(si, mean(auc_rep), std(auc_rep), std(auc_rep), 2, 0.2, [0 0 0]);
    X_lbl{si,1} = sprintf('%s-%s', res_info{3}, res_info{5}(5:7));
end
xlim([0 n_res+1]);
ylim([0.52 0.67]);
y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', X_lbl, 'XTickLabelRotation', 20, ...
	'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

%% Saving
% output_name = sprintf('./Plots/S02_SingleCls_PerformanceComparison.pdf');
% set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [10 3], 'PaperPosition', [0 0 10 3]);
% print('-dpdf', '-r300', output_name);

%% Compare p-values
close all
pval_mat = zeros(n_res);
for si=1:n_res
    for sj=1:n_res
        [h,p,ci,stats] = ttest2(AUC_cmb(:,si), AUC_cmb(:,sj));
        if stats.tstat>0
            pval_mat(si,sj) = -log10(p);
        else
            pval_mat(si,sj) = log10(p);
        end
    end
end
imagesc(pval_mat);
colormap(jet(n_res));
set(gca, 'XTick', 1:n_res, 'XTickLabel', X_lbl, 'XTickLabelRotation', 20, ...
    'YTick', 1:n_res, 'YTickLabel', X_lbl);
caxis([-5 5]);


