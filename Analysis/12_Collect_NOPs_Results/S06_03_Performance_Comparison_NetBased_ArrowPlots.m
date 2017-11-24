clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    {'MRK_CVT01_Lasso_Random-G00500_MSN-500.mat' 'MRK_CVT01_GLasso_Random-G00500_MSN-500.mat'}
    {'MRK_CVT01_Lasso_I2D-G00500_MSN-500.mat' 'MRK_CVT01_GLasso_I2D-G00500_MSN-500.mat'}
    {'MRK_CVT01_Lasso_MSigDB-G00500_MSN-500.mat' 'MRK_CVT01_GLasso_MSigDB-G00500_MSN-500.mat'}
    {'MRK_CVT01_Lasso_HPRD-G00500_MSN-500.mat' 'MRK_CVT01_GLasso_HPRD-G00500_MSN-500.mat'}
    {'MRK_CVT01_Lasso_KEGG-G00500_MSN-500.mat' 'MRK_CVT01_GLasso_KEGG-G00500_MSN-500.mat'}
    {'MRK_CVT01_Lasso_STRING-G00500_MSN-500.mat' 'MRK_CVT01_GLasso_STRING-G00500_MSN-500.mat'}
    {'MRK_CVT01_Lasso_HBGland-G11748_MSN-500.mat' 'MRK_CVT01_GLasso_HBGland-G11748_MSN-500.mat'}
    {'MRK_CVT01_Lasso_HBEpith-G11748_MSN-500.mat' 'MRK_CVT01_GLasso_HBEpith-G11748_MSN-500.mat'}
    {'MRK_CVT01_Lasso_ACr-G00500_MSN-500.mat' 'MRK_CVT01_GLasso_ACr-G00500_MSN-500.mat'}
    {'MRK_CVT01_Lasso_AvgSynACr-G00500_MSN-500.mat' 'MRK_CVT01_GLasso_AvgSynACr-G00500_MSN-500.mat'}
    };
n_res = numel(res_lst);
indData_path = 'MRK_CVT01_LExAG_None-G11748_MSN-500.mat';

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(n_res);
Method_lst = {};
xlim([0.595 0.675]);
ylim([0.015 0.045]);
for si=1:n_res
    [auc1_mat, ~       ] = LoadData([res_path res_lst{si}{1}]);
    [auc2_mat, res_info] = LoadData([res_path res_lst{si}{2}]);
    Method_lst{si, 1} = res_info{5};
    
    [bar_h(si,:), avg1_rep, std1_fld] = DrawHorzErrBar(auc1_mat, 0.001, 'Color', clr_map(si,:), 'LineWidth', 4);
    [bar_h(si,:), avg2_rep, std2_fld] = DrawHorzErrBar(auc2_mat, 0.001, 'Color', clr_map(si,:), 'LineWidth', 4);
    plot([avg1_rep avg2_rep], [std1_fld std2_fld], ':', 'LineWidth', 4, 'Color', clr_map(si,:));
end

auc_mat = LoadData([res_path indData_path]);
bar_h(si+1,:) = DrawHorzErrBar(auc_mat, 0.001, 'Color', [0 0 0], 'LineWidth', 4);
Method_lst(end-1:end+1) = {'Correlation' 'SyNet' 'No network (all genes)'};

legend(bar_h(:,1), Method_lst, 'FontSize', 10, 'FontWeight', 'Bold', 'Location', 'SouthEastOutside');

set(gca, 'XLim', [0.595 0.675], 'YLim', [0.015 0.043]);
set(gca, 'FontWeight', 'Bold', 'FontSize', 12);
xlabel('Average AUC (across repeats)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Standard deviation (across folds)', 'FontSize', 12, 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_NetBased_ArrowBased.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [12 4], 'PaperPosition', [0 0 12 4]);
print('-dpdf', '-r300', output_name);

%% Function %%%%%%%%%%%%%%%%%%%%
function [auc_mat, res_info] = LoadData(res_name)
fprintf('Reading [%s]\n', res_name);
res_data = load(res_name);
if any(isnan(res_data.AUC_mat(:))), error(); end
auc_mat = res_data.AUC_mat;
res_info = regexp(res_name, '[_-]', 'split');
end

function [bar_h, avg_rep, std_fld] = DrawHorzErrBar(auc_mat, wiskerwidth, varargin)
linewidth = 2;
color = [0 0 0];

auc_rep = mean(auc_mat,1);
std_rep = std(auc_rep, 0, 2);
avg_rep = mean(auc_rep);
std_fld = mean(std(auc_mat, 0, 2));

bar_h(1,1) = plot(avg_rep+[-std_rep +std_rep], std_fld([1 1]), 'linewidth', linewidth, 'color', color, varargin{:});
bar_h(1,2) = plot(avg_rep+[+std_rep +std_rep], std_fld+[-wiskerwidth wiskerwidth], 'linewidth', linewidth, 'color', color, varargin{:});
bar_h(1,3) = plot(avg_rep+[-std_rep -std_rep], std_fld+[-wiskerwidth wiskerwidth], 'linewidth', linewidth, 'color', color, varargin{:});
end