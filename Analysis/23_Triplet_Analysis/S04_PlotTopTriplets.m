clc;
clear;

%% Initialization
addpath('../_Utilities/', '-end');
Triplet_AUC = zeros(0, 11);

%% Get top results
res_ptr = sprintf('./TRC_Files/SyNet/TC_SyNet_*.mat');
res_lst = dir(res_ptr);
for ri=1:numel(res_lst)
    fprintf('Loading [%s] ...\n', res_lst(ri).name);
    res_name = [res_lst(ri).folder '/' res_lst(ri).name];
    res_info = load(res_name);
    Triplet_AUC = [Triplet_AUC; res_info.Triplet_AUC(1:1000,:)];
end

%% Sorting
[~, sid] = sort(Triplet_AUC(:,11), 'Descend');
Triplet_AUC = Triplet_AUC(sid,:);

%% Plotting
close all
clr_map = [
    0.7 0.7 0.7
    0.4 0.4 0.8
    0.9 0.3 0.2
    ];
figure('Position', [100 100 1400 600]);
hold on
n_top = 100;
for bi=4:6
    BoxPlotEx(Triplet_AUC(1:n_top, bi), 'Positions', bi/4-0.25, 'Color', clr_map(1,:), 'Width', 0.2);
end
for bi=7:9
    BoxPlotEx(Triplet_AUC(1:n_top, bi), 'Positions', bi/4, 'Color', clr_map(2,:), 'Width', 0.2);
end
BoxPlotEx(Triplet_AUC(1:n_top, 10), 'Positions', 3, 'Color', clr_map(3,:), 'Width', 0.2);
set(gca, 'XTick', 1:3, 'XTickLabel', {'Genes', 'Pairs', 'Triplets'}, 'FontWeight', 'Bold', ...
    'FontSize', 14);
xlim([0.5 3.5]);
ylim([0.53 0.65]);
title('Performance comparison of single/pair/triple wise genes');

%% Saving
output_name = sprintf('./Plots/S04_PerformanceOfTriplets_nTop%d.pdf', n_top);
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [12 5], 'PaperPosition', [0 0 12 5]);
print('-dpdf', '-r300', output_name);
