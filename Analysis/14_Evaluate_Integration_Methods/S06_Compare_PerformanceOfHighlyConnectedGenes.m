clc;
clear;

%% Initialization
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

%% Collection of data
IND_auc = zeros(n_gene, n_rep, n_study);
for ci=1:n_study
    for ri=1:n_rep
        res_name = sprintf('ResIND_SyNet-SyNet_CVT50_Si%02d-Ri%03d.mat', ci, ri);
        fprintf('Result: %s\n', res_name);
        res_data = load(['./Result_Files/' res_name]);
        IND_auc(:, ri, ci) = res_data.Gene_TrAUC;
    end
end
IND_auc = mean(mean(IND_auc, 3), 2);

%% Compute frequency of genes
str_freq = histcounts(horzcat(STR_nei{:}), 1:n_gene+1);
shf_freq = histcounts(horzcat(SHF_nei{:}), 1:n_gene+1);
bin_lst = floor(linspace(1,n_gene,7));

%% Plotting
close all
figure('Position', [100 100 1500 500]);
hold on
freq_set = floor(logspace(log10(1), log10(101), 7))';
freq_class = [freq_set(1:end-1) freq_set(2:end)-1];
freq_class(end) = inf;
n_class = size(freq_class, 1);
clr_map = [
    0.0 0.0 1.0;
    0.7 0.7 0.7;
    ];
x_label = cell(n_class, 1);
% IND_auc(IND_auc<0.55) = nan;
for ci=1:n_class
    has_ol = shf_freq>=freq_class(ci,1) & shf_freq<=freq_class(ci,2);
    box_h = boxplot(IND_auc(has_ol), 'Positions', ci-0.2, 'Color', clr_map(2,:));
    set(box_h, 'LineWidth', 2);
    
    has_ol = str_freq>=freq_class(ci,1) & str_freq<=freq_class(ci,2);
    box_h = boxplot(IND_auc(has_ol), 'Positions', ci+0.2, 'Color', clr_map(1,:));
    set(box_h, 'LineWidth', 2);
    
    text(ci, 0.50, sprintf('%d : %d\nn=%d', freq_class(ci,:), sum(has_ol)), 'HorizontalAlignment', 'Center');
end
xlim([0 n_class+1]);
ylim([0.50 0.60]);
set(gca, 'XTick', 1:n_class, 'XTickLabel', []);

%% Saving the plot
output_name = sprintf('./Plots/S06_PerformanceComparison_STRING_PerFreq.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [16 4], 'PaperPosition', [0 0 16 4]);
print('-dpdf', '-r300', output_name);

