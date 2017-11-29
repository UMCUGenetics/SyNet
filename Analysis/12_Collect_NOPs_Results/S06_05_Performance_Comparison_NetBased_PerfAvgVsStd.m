clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/SpiderPlot/');

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
figure('Position', [100 100 700 400]);
hold on
xlabel('Average AUC', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Standard deviation', 'FontSize', 12, 'FontWeight', 'Bold');
clr_map = jet(n_res)*0.95;
Method_lst = cell(n_res, 1);

for si=1:n_res
    [auc1_mat, ~       ] = LoadData([res_path res_lst{si}{1}]);
    [auc2_mat, res_info] = LoadData([res_path res_lst{si}{2}]);
    Method_lst{si, 1} = res_info{5};
    
    [avg1_tot, std1_rep, std1_study] = GetStdAvg(auc1_mat);
    [avg2_tot, std2_rep, std2_study] = GetStdAvg(auc2_mat);
    
    line_h(si,1) = plot([avg1_tot avg2_tot], [std1_rep std2_rep], 'Color', clr_map(si,:), 'LineWidth', 3);
    plot(avg2_tot, std2_rep, 'O', 'MarkerFaceColor', clr_map(si,:)*0.9, 'MarkerEdgeColor', clr_map(si,:)*0.8, 'MarkerSize', 10);
end
Method_lst(end-1:end+1) = {'Correlation' 'SyNet' ' '};
legend(line_h, Method_lst(1:end-1), 'FontWeight', 'Bold', 'Location', 'NorthEast');

%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_NetBased_PerfAvgVsStd.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'Landscape', 'PaperPositionMode','auto', 'PaperSize', [7 4], 'PaperPosition', [0 0 7 4]);
print('-dpdf', '-r300', output_name);

%% Function %%%%%%%%%%%%%%%%%%%%
function [auc_mat, res_info] = LoadData(res_name)
fprintf('Reading [%s]\n', res_name);
res_data = load(res_name);
if any(isnan(res_data.AUC_mat(:))), error(); end
auc_mat = res_data.AUC_mat;
res_info = regexp(res_name, '[_-]', 'split');
end

function [avg_total, std_repeat, std_study] = GetStdAvg(auc_mat)
avg_study = mean(auc_mat,1);
std_repeat = std(avg_study);
avg_total = mean(avg_study);
% std_study = mean(std(auc_mat, 0, 2));
% std_study = mean(std(auc_mat, 0, 2));
std_study = std(mean(auc_mat, 2));
end

