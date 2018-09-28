clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx');

%% Load nets
res_path = './SyNet_Overlap/';
net_lst = {
    'HumanInt' 'BioPlex'  'BioGRID' 'IntAct' 'STRING' 'HBOvary' 'HBBrain' 'HBKidney' 'HBLympNode' 'HBGland' ...
    };
n_net = numel(net_lst);
max_Y = 500;
Y_Scale = 'Linear';
y_lim = [0 max_Y];
SyNet_SIZE = 3544;
LIMIT_GENES = 0;
PLOT_SHUFFLE = 0;
if LIMIT_GENES
    limit_method = 'LimitedToRef';
else
    limit_method = 'All';
end
if PLOT_SHUFFLE
    data_src = '-SHFL';
else
    data_src = '';
end
MAX_N_PAIR_lst = [1000, 10000, 50000, 100000, 500000, 1000000];
n_pairs = numel(MAX_N_PAIR_lst);

%% Loop over max number of pairs
ovl_avg = zeros(n_net, n_pairs);
ovl_std = zeros(n_net, n_pairs);
ovl_clr = cell(n_net, n_pairs);
ovl_name = cell(n_net, n_pairs);
for ni=1:n_net
    for pi=1:n_pairs
        res_name = sprintf([res_path 'NetOV_SyNet_%s%s_MP%d_%s_SS%d.mat'], net_lst{ni}, data_src, MAX_N_PAIR_lst(pi), limit_method, SyNet_SIZE);
        %     if ~exist(res_name, 'file')
        %         fprintf('Warning: [%s] file is not found.\n', res_name);
        %         continue;
        %     end
        fprintf('Reading [%s]\n', res_name);
        res_data = load(res_name, 'OL_Freq');
        res_ovl = res_data.OL_Freq;
        % res_ovl(res_ovl>max_Y) = max_Y;
        
        ovl_avg(ni, pi) = mean(res_ovl);
        ovl_std(ni, pi) = std(res_ovl);
        
        res_info = regexp(res_name, '_', 'split');
        [met_clr, met_name] = getColor(res_info{4});
        ovl_name{ni, pi} = met_name;
        ovl_clr{ni, pi} = met_clr; 
    end
end

%% Plotting performance
close all
figure('Position', [100 100 700 400]);
hold on
met_hwnd = cell(n_net, 1);
for ni=1:n_net
    plot(1:n_pairs, ovl_avg(ni, :), 'Color', ovl_clr{ni, pi});
    for pi=1:n_pairs
        met_hwnd{ni, 1} = plot(pi, ovl_avg(ni, pi), 'Marker', 'O', 'MarkerSize', 7, 'MarkerFaceColor', ovl_clr{ni, pi}, 'MarkerEdgeColor', 'none');
        errorbarEx(pi, ovl_avg(ni, pi), ovl_std(ni, pi), ovl_std(ni, pi), 2, 0.1, ovl_clr{ni, pi});
    end
end
xlim([0 n_pairs+1]);
ylim(y_lim);
y_tick = linspace(0, max_Y, 11);
Y_lbl = arrayfun(@(y) sprintf('%0.0f', y), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_pairs, 'XTickLabel', [], 'XTickLabelRotation', 0, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'YScale', Y_Scale, 'YMinorTick','off', 'YMinorGrid','off', 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.15);
ylabel('# SyNet links', 'FontWeight', 'Bold');
title(sprintf('Overlap between networks limited on top links (%s)', data_src));

%% Saving
output_name = sprintf('./Plots/S07_OverlappingComparisonAcrossNumPairs_Src-SyNet_Use%s_%s_SS%d.pdf', data_src, Y_Scale, SyNet_SIZE);
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [12 4], 'PaperPosition', [0 0 12 4]);
print('-dpdf', '-r300', output_name);

% for ni=1:n_net
%     for pi=1:n_pairs
%             res_info = regexp(res_name, '_', 'split');
%             [met_clr, Method_lbl{bi,1}] = getColor(res_info{4});
%             ovl_clr{ni, pi} = met_clr;
%             box_h = BoxPlotEx(data_net, 'Positions', bi+offset, 'Color', met_clr, 'Symbol', '', 'Widths', 0.3);
%             set(box_h, 'LineWidth', 2);
%
%             %net_name = sprintf('%s\n#G%d, #L%d', Method_lbl{si}, numel(res_data.Net_GeneName), res_data.Net_nlnk);
%             net_name = sprintf('%s', Method_lbl{bi});
%             %net_name = sprintf('%s\n%s', Method_lbl{n_bar}, ref_lst{ri});
%             text(bi, y_lim(1)-1, net_name, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', ...
%                 'FontSize', 10, 'FontWeight', 'Bold', 'Rotation', 0);
%         end
%         bi = bi + 1;
%     end
%     bi = bi + 5;
% end
% xlim([0 bi+1]);
% ylim(y_lim);
% % y_tick = 2.^(0:11);
% y_tick = linspace(0, max_Y, 11);
% Y_lbl = arrayfun(@(y) sprintf('%0.0f', y), y_tick, 'UniformOutput', 0);
% set(gca, 'XTick', 1:bi, 'XTickLabel', [], 'XTickLabelRotation', 0, ...
%     'YTick', y_tick, 'YTickLabel', Y_lbl, 'YScale', Y_Scale, 'YMinorTick','off', 'YMinorGrid','off', 'FontWeight', 'Bold', 'FontSize', 10, ...
%     'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.15);
% ylabel('# SyNet links', 'FontWeight', 'Bold');
% title(sprintf('Overlap between networks limited on %d top links', MAX_N_PAIR))

% return
%% Saving
% output_name = sprintf('./Plots/S07_OverlappingComparisonAcrossNumPairs_%s_%s_SS%d.pdf', ref_lst{1}, Y_Scale, SAMPLE_SIZE);
% set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [12 4], 'PaperPosition', [0 0 12 4]);
% print('-dpdf', '-r300', output_name);



