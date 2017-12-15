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
%     {'NetOV_HBSpermatocyte-SHFL_NL10000.mat' 'NetOV_HBSpermatocyte_NL10000.mat'}
%     {'NetOV_HBSpermatogonium-SHFL_NL10000.mat' 'NetOV_HBSpermatogonium_NL10000.mat'}
%     {'NetOV_HBSpermatid-SHFL_NL10000.mat' 'NetOV_HBSpermatid_NL10000.mat'}
    {'NetOV_BioGRID-SHFL_NL10000.mat' 'NetOV_BioGRID_NL10000.mat'}
    {'NetOV_BioPlex-SHFL_NL10000.mat' 'NetOV_BioPlex_NL10000.mat'}
    {'NetOV_IntAct-SHFL_NL10000.mat' 'NetOV_IntAct_NL10000.mat'}
%     {'NetOV_HBNeuron-SHFL_NL10000.mat' 'NetOV_HBNeuron_NL10000.mat'}
    {'NetOV_STRING-SHFL_NL10000.mat' 'NetOV_STRING_NL10000.mat'}
%     {'NetOV_HBEpith-SHFL_NL10000.mat' 'NetOV_HBEpith_NL10000.mat'}
%     {'NetOV_HBRetina-SHFL_NL10000.mat' 'NetOV_HBRetina_NL10000.mat'}
%     {'NetOV_HBStomach-SHFL_NL10000.mat' 'NetOV_HBStomach_NL10000.mat'}
    {'NetOV_HBBrain-SHFL_NL10000.mat' 'NetOV_HBBrain_NL10000.mat'}
%     {'NetOV_HBNervs-SHFL_NL10000.mat' 'NetOV_HBNervs_NL10000.mat'}
    {'NetOV_HBLiver-SHFL_NL10000.mat' 'NetOV_HBLiver_NL10000.mat'}
%     {'NetOV_HBMuscle-SHFL_NL10000.mat' 'NetOV_HBMuscle_NL10000.mat'}
%     {'NetOV_HBHeart-SHFL_NL10000.mat' 'NetOV_HBHeart_NL10000.mat'}
%     {'NetOV_HBKidney-SHFL_NL10000.mat' 'NetOV_HBKidney_NL10000.mat'}
%     {'NetOV_HBOvary-SHFL_NL10000.mat' 'NetOV_HBOvary_NL10000.mat'}
%     {'NetOV_HBLung-SHFL_NL10000.mat' 'NetOV_HBLung_NL10000.mat'}
%     {'NetOV_HBLympNode-SHFL_NL10000.mat' 'NetOV_HBLympNode_NL10000.mat'}
%     {'NetOV_HBProstate-SHFL_NL10000.mat' 'NetOV_HBProstate_NL10000.mat'}
%     {'NetOV_HBEye-SHFL_NL10000.mat' 'NetOV_HBEye_NL10000.mat'}
%     {'NetOV_HBBone-SHFL_NL10000.mat' 'NetOV_HBBone_NL10000.mat'}
%     {'NetOV_HBColon-SHFL_NL10000.mat' 'NetOV_HBColon_NL10000.mat'}
%     {'NetOV_HBTestis-SHFL_NL10000.mat' 'NetOV_HBTestis_NL10000.mat'}
%     {'NetOV_HBBlood-SHFL_NL10000.mat' 'NetOV_HBBlood_NL10000.mat'}
%     {'NetOV_HBIntestine-SHFL_NL10000.mat' 'NetOV_HBIntestine_NL10000.mat'}
%     {'NetOV_HBEpidermis-SHFL_NL10000.mat' 'NetOV_HBEpidermis_NL10000.mat'}
    {'NetOV_HBGland-SHFL_NL10000.mat' 'NetOV_HBGland_NL10000.mat'}
%     {'NetOV_HBUterus-SHFL_NL10000.mat' 'NetOV_HBUterus_NL10000.mat'}
    % {'NetOV_AbsCorr-SHFL_NL10000.mat' 'NetOV_AbsCorr_NL10000.mat'}
};
n_res = numel(res_lst);
clr_map = lines(n_res);
max_Y = 2048;
Y_Scale = 'Log';

%% Plotting performance
close all
figure('Position', [100 100 700 400]);
hold on
X_lbl = {};
y_lim = [0.9 max_Y];
offset = 0.15;
for si=1:n_res
    res_name = [res_path res_lst{si}{1}];
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    data_rnd = res_data.OL_Freq;
    data_rnd(data_rnd>max_Y) = max_Y;
    
    box_h = BoxPlotEx(data_rnd, 'Positions', si-offset, 'Color', [0.7 0.7 0.7], 'Symbol', '', 'Widths', 0.3);
    set(box_h, 'LineWidth', 2);
    
    res_name = [res_path res_lst{si}{2}];
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    data_net = res_data.OL_Freq;
    data_net(data_net>max_Y) = max_Y;
    
    res_info = regexp(res_lst{si}{2}, '_', 'split');
    [met_clr, Method_lbl{si,1}] = getColor(res_info{2});
    box_h = BoxPlotEx(data_net, 'Positions', si+offset, 'Color', met_clr, 'Symbol', '', 'Widths', 0.3);
    set(box_h, 'LineWidth', 2);
    
    [~, pval] = ttest2(data_rnd, data_net);
    
    %net_name = sprintf('%s\n#G%d, #L%d', Method_lbl{si}, numel(res_data.Net_GeneName), res_data.Net_nlnk);
    net_name = sprintf('%s', Method_lbl{si});
    text(si, y_lim(1), net_name, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', ...
        'FontSize', 10, 'FontWeight', 'Bold');
end
% legend([left_h{1} right_h{1}], {'Random', 'STRING'});
xlim([0 n_res+1]);
ylim(y_lim);
y_tick = 2.^(0:11);
% y_tick = round(linspace(0, max_Y, 11));
Y_lbl = arrayfun(@(y) sprintf('%0.0f', y), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', 1:n_res, 'XTickLabel', [], 'XTickLabelRotation', 0, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'YScale', Y_Scale, 'YMinorTick','off', 'YMinorGrid','off', 'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.15);
ylabel('# SyNet links', 'FontWeight', 'Bold');

% return
%% Saving
output_name = sprintf('./Plots/S02_OverlappingComparison_%s.pdf', Y_Scale);
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [10 4], 'PaperPosition', [0 0 10 4]);
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