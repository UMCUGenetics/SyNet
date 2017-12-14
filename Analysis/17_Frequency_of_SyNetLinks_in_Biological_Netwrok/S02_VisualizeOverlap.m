clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../_Utilities/');

%% Load nets
res_path = './SyNet_Overlap/';
res_lst = {
    {'NetOV_BioPlex-SHFL_NL10000.mat' 'NetOV_BioPlex_NL10000.mat'}
    {'NetOV_BioGRID-SHFL_NL10000.mat' 'NetOV_BioGRID_NL10000.mat'}
    {'NetOV_IntAct-SHFL_NL10000.mat' 'NetOV_IntAct_NL10000.mat'}
    {'NetOV_STRING-SHFL_NL10000.mat' 'NetOV_STRING_NL10000.mat'}
    {'NetOV_HBBone-SHFL_NL10000.mat' 'NetOV_HBBone_NL10000.mat'}
    {'NetOV_HBBrain-SHFL_NL10000.mat' 'NetOV_HBBrain_NL10000.mat'}
    {'NetOV_HBEpith-SHFL_NL10000.mat' 'NetOV_HBEpith_NL10000.mat'}
    {'NetOV_HBLympNode-SHFL_NL10000.mat' 'NetOV_HBLympNode_NL10000.mat'}
    {'NetOV_HBGland-SHFL_NL10000.mat' 'NetOV_HBGland_NL10000.mat'}
%     {'NetOV_AbsCorr-SHFL_NL10000.mat' 'NetOV_AbsCorr_NL10000.mat'}
};
n_res = numel(res_lst);
clr_map = lines(n_res);
max_Y = 1500;
edge_lst = floor(linspace(0, max_Y, 10));

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
X_lbl = {};
y_lim = [0.9 max(edge_lst)];
offset = 0.15;
for si=1:n_res
    res_name = [res_path res_lst{si}{1}];
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    data_rnd = res_data.OL_Freq;
    data_rnd(data_rnd>max_Y) = max_Y;
    
    %left_h = ViolinEx(si-offset, edge_lst, data_rnd, struct('ShrinkFactor', 0.5, 'BarColor', [0.7 0.7 0.7], 'Reverse', 1));
    %set(left_h, 'FaceAlpha', 0.8);
    box_h = BoxPlotEx(data_rnd, 'Positions', si-offset, 'Color', [0.7 0.7 0.7], 'Symbol', '', 'Widths', 0.25);
    set(box_h, 'LineWidth', 2);
    
    res_name = [res_path res_lst{si}{2}];
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    data_net = res_data.OL_Freq;
    data_net(data_net>max_Y) = max_Y;
    
    %right_h = ViolinEx(si+offset, edge_lst, data_net, struct('ShrinkFactor', 0.5, 'BarColor', clr_map(si,:), 'Reverse', 0));
    %right_h = distributionPlot(data_net, 'color', clr_map(si,:), 'xValues', si+offset+0.4, 'showMM', 0, 'distWidth', 0.7, 'histOri', 'right', 'divFactor', 1);
    %set(right_h, 'FaceAlpha', 0.8);
    res_info = regexp(res_lst{si}{2}, '_', 'split');
    [met_clr, Method_lbl{si,1}] = getColor(res_info{2});
    box_h = BoxPlotEx(data_net, 'Positions', si+offset, 'Color', met_clr, 'Symbol', '', 'Widths', 0.25);
    set(box_h, 'LineWidth', 2);
    
    [~, pval] = ttest2(data_rnd, data_net);
    
    net_name = sprintf('%s\n#G%d, #L%d', Method_lbl{si}, numel(res_data.Net_GeneName), res_data.Net_nlnk);
    text(si, y_lim(1), net_name, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', ...
        'FontSize', 12, 'FontWeight', 'Bold');
end
% legend([left_h{1} right_h{1}], {'Random', 'STRING'});
xlim([0 n_res+1]);
ylim(y_lim);
y_tick = 2.^(0:10); %[0 ceil(logspace(log10(1), log10(120), 8))]; %get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f', y), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', [], 'XTickLabelRotation', 0, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'YScale', 'Log', 'YMinorTick','off', 'YMinorGrid','off', 'FontWeight', 'Bold', 'FontSize', 12, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.15);
ylabel('# SyNet links', 'FontWeight', 'Bold');

return
%% Saving
output_name = sprintf('./Plots/S02_OverlappingComparison.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [16 4], 'PaperPosition', [0 0 16 4]);
print('-dpdf', '-r300', output_name);


function [patch_h, bin_freq] = ViolinEx(x_Pos, edge_lst, dist, opt)
if ~exist('opt', 'var'), opt = struct; end
if ~isfield(opt, 'BarColor')
    opt.BarColor = [0.1 0.1 1];
end
n_bin = numel(edge_lst)-1;
hold on

%% Make PDF
bin_freq = zeros(1, n_bin);
for bi=1:n_bin
    if bi==n_bin
        has_ol = edge_lst(bi)<=dist & dist <= edge_lst(bi+1);
    else
        has_ol = edge_lst(bi)<=dist & dist <  edge_lst(bi+1);
    end
    bin_freq(bi) = sum(has_ol);
end

%% Normalize and adjustments
if ~isfield(opt, 'Normalize') || opt.Normalize==1
    bin_nrm = bin_freq/max(bin_freq);
else
    bin_nrm = bin_freq;
end
if isfield(opt, 'ShrinkFactor')
    bin_nrm = bin_nrm * opt.ShrinkFactor;
end
if isfield(opt, 'Reverse') && opt.Reverse==1
    bin_nrm = - bin_nrm;
end

%% Plot bars
for bi=1:n_bin
    patch_h(bi,1) = patch([x_Pos x_Pos+bin_nrm([bi bi]) x_Pos], edge_lst([bi bi bi+1 bi+1]), 'b', 'EdgeColor', 'None', 'FaceColor', opt.BarColor);
end
end