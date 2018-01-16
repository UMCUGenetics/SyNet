clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
Ratio_lst = 0.05:0.1:1;
n_ratio = numel(Ratio_lst);
MAX_N_PAIR = 50000;
MAX_SyNet_Pairs = 03544;
MAX_SyNet_Genes = 300;
Net_lst = {'HBGland','HBLympNode','HBOvary','HBKidney','HBBrain','STRING','IntAct','BioGRID','BioPlex','HumanInt'}; % ,'HBGland-SHFL'
n_net = numel(Net_lst); 
Mrk_lst = {'O' 'x' '>' '<' '^' 'v' '+' 's' 'd' 'p' '-'};

% Ref_Name = 'AvgSyn';
Ref_Name = 'SyNet'; 
Limit_method = 'All';
% Limit_method = 'LimitedToRef';
% Overlap_Check = 'Gene';
Overlap_Check = 'Lnk';

%% Collecting data
pr_mat = zeros(n_net, n_ratio);
rc_mat = zeros(n_net, n_ratio);
for ni=1:n_net
    for ri=1:n_ratio
        res_name = sprintf('./SyNet_OverlapPerThresh/%sOV-PerThresh_%s_%s_MP%d_%s_RV%0.2f.mat', Overlap_Check, Ref_Name, Net_lst{ni}, MAX_N_PAIR, Limit_method, Ratio_lst(ri));
        fprintf('Loading [%s]\n', res_name);
        
        if strcmp(Overlap_Check, 'Gene')
            res_data = load(res_name, 'Real_Freq', 'n_SelGene');
            pr_mat(ni,ri) = res_data.Real_Freq / res_data.n_SelGene;
            rc_mat(ni,ri) = res_data.Real_Freq / MAX_SyNet_Genes;
        else
            res_data = load(res_name, 'Real_Freq', 'n_TopPair');
            pr_mat(ni,ri) = res_data.Real_Freq / res_data.n_TopPair;
            rc_mat(ni,ri) = res_data.Real_Freq / MAX_SyNet_Pairs;
        end
    end
end


%% Plotting
figure('Position', [100 100 1400 500]);
hold on
for ni=1:n_net
    [clr_map, Method_lst{ni,1}] = getColor(Net_lst{ni});
    plt_h = plot(rc_mat(ni,:), pr_mat(ni,:), [Mrk_lst{ni} '-'], 'Color', clr_map, 'MarkerFaceColor', 'None', ...
        'MarkerSize', 2, 'MarkerEdgeColor', clr_map, 'LineWidth', 2);
    plt_h.Annotation.LegendInformation.IconDisplayStyle = 'off';
    
    f1_scr = 2 * (pr_mat(ni,:) .* rc_mat(ni,:)) ./ (pr_mat(ni,:) + rc_mat(ni,:));
    [~, best_ind] = max(f1_scr);
    mrk_h = plot(rc_mat(ni,best_ind), pr_mat(ni,best_ind), Mrk_lst{ni}, 'MarkerFaceColor', clr_map * 0.8, ...
        'MarkerSize', 8, 'MarkerEdgeColor', clr_map * 0.8, 'LineWidth', 2);
    %mrk_h.Annotation.LegendInformation.IconDisplayStyle = 'off';
end
% xlim([0 1]);
% ylim([0 1]);
set(gca, 'FontWeight', 'Bold', 'FontSize', 12);
xlabel('Recall (sensitivity)');
ylabel('Precision');
title(sprintf('Overlap of [%s] and biological networks [Over %s, Method = %s]', Ref_Name, Overlap_Check, Limit_method));
if strcmp(Overlap_Check, 'Gene')
    legend(Method_lst);
end

% return
%% Save
output_name = sprintf('./Plots/S06_NetOL_RF-%s_OV-%s_LM-%s_MP%06d_PresRecall.pdf', Ref_Name, Overlap_Check, Limit_method, MAX_N_PAIR);
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [7 5], 'PaperPosition', [0 0 7 5]);
print('-dpdf', '-r300', output_name);


