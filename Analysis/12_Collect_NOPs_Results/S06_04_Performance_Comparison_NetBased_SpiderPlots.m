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
figure('Position', [100 100 1500 400]);
subplot(1, 2, 2);
hold on
xlabel('Standard deviation (across repeats)', 'FontSize', 12, 'FontWeight', 'Bold');
ylabel('Standard deviation (across study)', 'FontSize', 12, 'FontWeight', 'Bold');
clr_map = jet(n_res)*0.95;
Method_lst = cell(n_res, 1);
% set(gca, 'XLim', [0.595 0.675], 'YLim', [0.015 0.043], 'FontWeight', 'Bold', 'FontSize', 12);
Spider_Data = zeros(2, n_res);
for si=1:n_res
    [auc1_mat, ~       ] = LoadData([res_path res_lst{si}{1}]);
    [auc2_mat, res_info] = LoadData([res_path res_lst{si}{2}]);
    Method_lst{si, 1} = res_info{5};
    
    [avg1_rep, std1_rep, std1_study] = GetStdAvg(auc1_mat);
    [avg2_rep, std2_rep, std2_study] = GetStdAvg(auc2_mat);
    Spider_Data(:, si) = [avg1_rep; avg2_rep];
    
    line_h(si,1) = plot([std1_rep std2_rep], [std1_study std2_study], 'Color', clr_map(si,:), 'LineWidth', 3);
    plot(std2_rep, std2_study, 'O', 'MarkerFaceColor', clr_map(si,:)*0.9, 'MarkerEdgeColor', clr_map(si,:)*0.8, 'MarkerSize', 10);
end
Method_lst(end-1:end+1) = {'Correlation' 'SyNet' ' '};
legend(line_h, Method_lst(1:end-1), 'FontWeight', 'Bold', 'Location', 'SouthEast');

sp_h = subplot(1, 2, 1);
Spider_Data(:, end+1) = 0.61;
Spider_Data(3, :) = ones(1, n_res+1)*0.61;
Spider_Data(4, :) = ones(1, n_res+1)*0.68;
spider_plot(Spider_Data*100, Method_lst, 7, 1,...
    'Marker', 'o',...
    'LineStyle', '-',...
    'LineWidth', 2,...
    'MarkerSize', 5);
view(-122.5, 90);
line_h = findobj(sp_h, 'Type', 'Line');
set(line_h(1), 'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.6 0.6 0.6], 'MarkerEdgeColor', 'none');
% set(line_h(1), 'Color', [0.9 0.9 0.9 0.0], 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none');
set(line_h(2), 'Color', [0.9 0.9 0.9 0.0], 'MarkerFaceColor', 'none', 'MarkerEdgeColor', 'none');
text_h = findobj(sp_h, 'Type', 'text');
for i=2:19
    txt_pos = get(text_h(i), 'Position');
    set(text_h(i), 'FontWeight', 'Bold', 'LineStyle', 'none', 'Background', 'none', 'Position', txt_pos*1.1, 'HorizontalAlignment', 'Center');
end
delete(text_h([1 20:end]));
% auc_mat = LoadData([res_path indData_path]);
% bar_h(si+1,:) = DrawHorzErrBar(auc_mat, 0.001, 'Color', [0 0 0], 'LineWidth', 4);

legend({'Individual genes' 'Network based'}, 'FontWeight', 'Bold', 'Location', 'SouthWest');

%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_NetBased_SpiderPlot.pdf');
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

function [avg_rep, std_rep, std_study] = GetStdAvg(auc_mat)
avg_study = mean(auc_mat,1);
std_rep = std(avg_study, 0, 2);
avg_rep = mean(avg_study);
std_study = mean(std(auc_mat, 0, 2));
end