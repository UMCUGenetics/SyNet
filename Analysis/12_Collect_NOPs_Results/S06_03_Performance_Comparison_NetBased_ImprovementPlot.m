clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/ColumnLegend/');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    {'MRK_CVT01_LExAG_None-G11748_MSN-500.mat' 'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HumanInt-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_HumanInt-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_BioPlex-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_BioPlex-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_BioGRID-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_BioGRID-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_IntAct-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_IntAct-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_STRING-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_STRING-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBOvary-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_HBOvary-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBBrain-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_HBBrain-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBKidney-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_HBKidney-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBGland-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_HBGland-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBLympNode-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_HBLympNode-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_ACr-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_ACr-P50000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_AvgSynACr-P50000_MSN-500.mat' 'MRK_CVT01_NetGL_AvgSynACr-P50000_MSN-500.mat'}
    };
n_res = numel(res_lst);
mrk_lst = {'x','d', 'o', '*', '<', 'v', '+', 'x', 's', '^', '>', '<', 'p'};
n_mrk = numel(mrk_lst);

%% Plotting performance
close all
figure('Position', [100 100 700 400]);
hold on
Method_lst = cell(n_res, 1);
plot([-5 5], [0  0], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
plot([ 0 0], [-5 5], ':', 'Color', [0.5 0.5 0.5], 'LineWidth', 1);
for si=1:n_res
    [auc1_mat, ~       ] = LoadData([res_path res_lst{si}{1}]);
    [auc2_mat, res_info] = LoadData([res_path res_lst{si}{2}]);
    [met_clr, Method_lst{si,1}] = getColor(res_info{5});
    
    [avg1_tot, std1_rep, std1_study] = GetStdAvg(auc1_mat);
    [avg2_tot, std2_rep, std2_study] = GetStdAvg(auc2_mat);
    if si==1
        ref_avg = avg2_tot;
        ref_std = std2_rep;
    end
    
    x_avgi = avg2_tot - ref_avg;
    y_stdi = ref_std - std2_rep;
    if ismember(mrk_lst{si}, {'x' '+'})
        mrk_size = 200;
    else
        mrk_size = 80;
    end
    mrk_h(si,1) = scatter(x_avgi, y_stdi, mrk_size, 'Filled', mrk_lst{si}, 'MarkerFaceColor', met_clr, 'MarkerFaceAlpha', 1.0, ...
        'MarkerEdgeColor', met_clr, 'MarkerEdgeAlpha', 1.0, 'LineWidth', 2);
end
xlim([-3.5 3]);
ylim([-0.2 0.8]);

legend(mrk_h, Method_lst, 'FontWeight', 'Bold', 'FontSize', 10, 'Location', 'SouthEast', 'Orientation', 'Vertical');
% columnlegend(mrk_h, 2, Method_lst);
set(gca, 'FontWeight', 'Bold', 'FontSize', 8);
xlabel('AUC Increase', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Stability Improvement', 'FontSize', 12, 'FontWeight', 'Bold');

% return
%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_03_NetBased_Improvement.pdf');
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

