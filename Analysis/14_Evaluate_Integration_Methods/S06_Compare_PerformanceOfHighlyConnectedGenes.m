clc;
clear;

%% Initialization
addpath('../_Utilities/');
n_study = 14;
n_rep = 5;

%% Infer shuffled genes
net_STRING = load('./NetNei_Files/NetNei_STRING_NN20.mat');
net_Shuffl = load('./NetNei_Files/NetNei_Shf-STRING_NN20.mat');
if ~isequal(net_STRING.Neig_cell, net_Shuffl.Neig_cell), error(); end
n_gene = numel(net_STRING.Gene_Name);
rnd_ID = net_Shuffl.rnd_ID;

%% Get gene degree
STR_nei = net_STRING.Neig_cell;
SHF_nei = cell(n_gene, 1);
for gi=1:n_gene
    SHF_nei{gi}= rnd_ID(net_STRING.Neig_cell{gi});
end
str_degree = cellfun('length', STR_nei);
shf_degree = cellfun('length', SHF_nei);
if ~isequal(str_degree, shf_degree), exist(); end

%% ttest
GEPath = getPath('SyNet');
load(GEPath, 'Gene_Expression', 'Patient_Label');
zData = zscore(Gene_Expression);
IND_pvl = -log10(ttest2Ex(zData, Patient_Label));

%% Compute frequency of genes
str_freq = histcounts(horzcat(STR_nei{:}), 1:n_gene+1);
shf_freq = histcounts(horzcat(SHF_nei{:}), 1:n_gene+1);

%% Plotting
close all
figure('Position', [100 100 1500 500]);
hold on
freq_set = unique(floor(logspace(log10(1), log10(100), 7)))';
freq_class = [freq_set(1:end-1) freq_set(2:end)-1];
freq_class(end) = inf;
n_class = size(freq_class, 1);
clr_map = [
    0.7 0.7 0.7;
    0.1 0.1 1.0;
    ];
for ci=1:n_class
    str_ol = shf_freq>=freq_class(ci,1) & shf_freq<=freq_class(ci,2);
    BoxPlotEx(IND_pvl(str_ol), 'Positions', ci-0.25, 'Color', clr_map(1,:), 'Symbol', '');
    
    shf_ol = str_freq>=freq_class(ci,1) & str_freq<=freq_class(ci,2);
    BoxPlotEx(IND_pvl(shf_ol), 'Positions', ci+0.25, 'Color', clr_map(2,:), 'Symbol', '');
    
    text(ci, 0, sprintf('Interval %d : %d\n#gene=%d', freq_class(ci,:), sum(shf_ol)), 'FontWeight', 'Bold', 'FontSize', 10, ...
    'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top');

    [~, pval] = ttest2(IND_pvl(shf_ol), IND_pvl(str_ol), 'Tail', 'Right');
    text(ci, 3, sprintf('p = %0.0e', pval), 'FontWeight', 'Bold', 'FontSize', 10, 'Color', 'r', ...
        'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Bottom');
end
title('Performance of genes according to their membership frequency in sub-networks', 'FontSize', 14);
line_h(1) = plot([-1 -1], 'Color', clr_map(1,:), 'LineWidth', 7);
line_h(2) = plot([-1 -1], 'Color', clr_map(2,:), 'LineWidth', 7);
legend(line_h, {'Shuffled STRING' 'STRING'}, 'FontWeight', 'Bold', 'Location', 'NorthWest');
ylabel('-Log10(p-value)', 'FontWeight', 'Bold');
xlim([0.5 n_class+0.5]);
ylim([0 9])
set(gca, 'XTick', 1:n_class, 'XTickLabel', []);

%% Saving the plot
output_name = sprintf('./Plots/S06_PerformanceComparison_STRING_PerFreq.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [16 4], 'PaperPosition', [0 0 16 4]);
print('-dpdf', '-r300', output_name);

