clc;
clear;

%% Initialization
addpath('../_Utilities/');
top_src = '../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat';
n_top = 10000;

%% Load data
top_data = load(top_src);
top_arr = top_data.PP_Info(1:n_top, [5 7]);
% 5: Maximum AUC
% 6: Combined AUC
% 7: Synergy
% 8: Mean AUC
% 9: Corr

%% Generate classes
cls_lst = prctile(top_arr(:,1), 0:20:100)';
cls_arr = [cls_lst(1:end-1) cls_lst(2:end)-1e-5];
% cls_arr([1 end]) = [-inf inf];
n_class = size(cls_arr, 1);

%% Plotting
close all
figure('Position', [100 100 1500 500]);
hold on
clr_map = winter(n_class);
y_lim = [0.94 1.07];
for ci=1:n_class
    has_ol = top_arr(:,1)>=cls_arr(ci,1) & top_arr(:,1)<=cls_arr(ci,2);
    if sum(has_ol)==0, continue; end
    BoxPlotEx(top_arr(has_ol,2), 'Widths', 0.7, 'Positions', ci, 'Color', clr_map(ci,:), 'Symbol', ''); % 
    
    text(ci, y_lim(1), sprintf('[%0.1f : %0.1f)', cls_arr(ci,:)*100), 'FontSize', 12, 'FontWeight', 'Bold', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');
end

title('Synergy gained according to maximum individual AUC of each pair', 'FontSize', 14);
% line_h(1) = plot([-1 -1], 'Color', clr_map(1,:), 'LineWidth', 7);
% line_h(2) = plot([-1 -1], 'Color', clr_map(2,:), 'LineWidth', 7);
% legend(line_h, {'Shuffled STRING' 'STRING'}, 'FontWeight', 'Bold', 'Location', 'NorthWest');
ylabel('Synergy', 'FontWeight', 'Bold');
xlim([0 n_class+1]);
ylim(y_lim);
set(gca, 'XTick', 1:n_class, 'XTickLabel', []);

%% Saving the plot
output_name = sprintf('./Plots/S07_Synergy_vs_MaxIndAUC.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [10 4], 'PaperPosition', [0 0 10 4]);
print('-dpdf', '-r300', output_name);

