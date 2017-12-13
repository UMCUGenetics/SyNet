clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/SpiderPlot/');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    {'MRK_CVT01_NetLasso_BioPlex-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_BioPlex-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_BioGRID-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_BioGRID-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_IntAct-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_IntAct-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_STRING-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_STRING-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBGland-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBGland-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBLympNode-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBLympNode-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBEpith-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBEpith-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_ACr-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_ACr-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_AvgSynACr-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_AvgSynACr-P25000_MSN-500.mat'}
    };
n_res = numel(res_lst);

%% Plotting performance
close all
figure('Position', [100 100 700 400]);
hold on
xlabel('Average AUC', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Standard deviation', 'FontSize', 12, 'FontWeight', 'Bold');
Method_lst = cell(n_res, 1);

for si=1:n_res
    [auc1_mat, ~       ] = LoadData([res_path res_lst{si}{1}]);
    [auc2_mat, res_info] = LoadData([res_path res_lst{si}{2}]);
    [met_clr, Method_lst{si,1}] = getColor(res_info{5});
    
    [avg1_tot, std1_rep, std1_study] = GetStdAvg(auc1_mat);
    [avg2_tot, std2_rep, std2_study] = GetStdAvg(auc2_mat);
    
    line_h(si,1) = plot([avg1_tot avg2_tot], [std1_rep std2_rep], 'Color', met_clr, 'LineWidth', 3);
    plot(avg1_tot, std1_rep, 'S', 'MarkerFaceColor',     [1 1 1], 'MarkerEdgeColor', met_clr*1.0, 'MarkerSize', 10, 'LineWidth', 2);
    plot(avg2_tot, std2_rep, 'O', 'MarkerFaceColor', met_clr*1.0, 'MarkerEdgeColor', met_clr*1.0, 'MarkerSize', 10);
end
legend(line_h, Method_lst, 'FontWeight', 'Bold', 'Location', 'NorthEast');

%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_NetBased_PerfAvgVsStd.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'Landscape', 'PaperPositionMode','auto', 'PaperSize', [5 4], 'PaperPosition', [0 0 5 4]);
print('-dpdf', '-r300', output_name);

%% Function %%%%%%%%%%%%%%%%%%%%
function [auc_mat, res_info] = LoadData(res_name)
fprintf('Reading [%s]\n', res_name);
res_data = load(res_name);
if any(isnan(res_data.AUC_mat(:))), error(); end
auc_mat = res_data.AUC_mat * 100;
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

