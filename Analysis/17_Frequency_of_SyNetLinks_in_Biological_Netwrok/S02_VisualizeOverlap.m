clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../_Utilities/');

%% Load nets
res_path = './SyNet_Overlap/';
net_lst = {'BioGRID' 'BioPlex' 'IntAct' 'STRING' 'HBBrain' 'HBLiver' 'HBEpith' 'HBLympNode' 'HBGland'};
n_net = numel(net_lst);
ref_lst = {'SyNet' 'AvgSynACr' 'ACr' 'Syn' 'Avg' 'AvgSyn'};
n_ref = numel(ref_lst);
max_Y = 2048;
Y_Scale = 'Log';
y_lim = [0.9 max_Y];
offset = 0.15;

%% Plotting performance
close all
figure('Position', [100 100 700 400]);
hold on
n_bar = 0;
for ni=1:n_net
    for ri=1:n_ref
        res_name = sprintf([res_path 'NetOV_%s_%s-SHFL_NL10000.mat'], ref_lst{ri}, net_lst{ni});
        if ~exist(res_name, 'file')
            fprintf('Warning: [%s] file is not found.\n', res_name);
            continue;
        end
        n_bar = n_bar + 1;
        fprintf('Reading [%s]\n', res_name);
        res_data = load(res_name, 'OL_Freq');
        data_rnd = res_data.OL_Freq;
        data_rnd(data_rnd>max_Y) = max_Y;
        
        box_h = BoxPlotEx(data_rnd, 'Positions', n_bar-offset, 'Color', [0.7 0.7 0.7], 'Symbol', '', 'Widths', 0.3);
        set(box_h, 'LineWidth', 2);
        
        res_name = sprintf([res_path 'NetOV_%s_%s_NL10000.mat'], ref_lst{ri}, net_lst{ni});
        fprintf('Reading [%s]\n', res_name);
        res_data = load(res_name, 'OL_Freq');
        data_net = res_data.OL_Freq;
        data_net(data_net>max_Y) = max_Y;
        
        res_info = regexp(res_name, '_', 'split');
        [met_clr, Method_lbl{n_bar,1}] = getColor(res_info{4});
        box_h = BoxPlotEx(data_net, 'Positions', n_bar+offset, 'Color', met_clr, 'Symbol', '', 'Widths', 0.3);
        set(box_h, 'LineWidth', 2);
        
        %net_name = sprintf('%s\n#G%d, #L%d', Method_lbl{si}, numel(res_data.Net_GeneName), res_data.Net_nlnk);
        %net_name = sprintf('%s', Method_lbl{n_bar});
        net_name = sprintf('%s\n%s', Method_lbl{n_bar}, ref_lst{ri});
        text(n_bar, y_lim(1), net_name, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', ...
            'FontSize', 10, 'FontWeight', 'Bold');
    end
end

xlim([0 n_bar+1]);
ylim(y_lim);
y_tick = 2.^(0:11);
Y_lbl = arrayfun(@(y) sprintf('%0.0f', y), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_bar, 'XTickLabel', [], 'XTickLabelRotation', 0, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'YScale', Y_Scale, 'YMinorTick','off', 'YMinorGrid','off', 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.15);
ylabel('# SyNet links', 'FontWeight', 'Bold');

return
%% Saving
output_name = sprintf('./Plots/S02_OverlappingComparison_%s.pdf', Y_Scale);
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [10 4], 'PaperPosition', [0 0 10 4]);
print('-dpdf', '-r300', output_name);



