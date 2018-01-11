clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
Ratio_lst = 0.05:0.05:1;
n_ratio = numel(Ratio_lst);
MAX_N_PAIR = 50000;
Ref_Name = 'SyNet'; %'AvgSyn'
Net_lst = {'HumanInt' 'BioPlex' 'BioGRID' 'IntAct' 'STRING' 'HBBrain' 'HBKidney' 'HBOvary' 'HBLympNode' 'HBGland'};
n_net = numel(Net_lst);
Mrk_lst = {'-' 'O' 'x' '>' '<' '^' 'v' '+' 's' 'd' 'p'};
Limit_method = 'All';
% Limit_method = 'LimitedToRef';
Overlap_Check = 'Gene';
% Overlap_Check = 'Lnk';

%% Collecting data
zscr_Mat = zeros(n_net, n_ratio);
param_str = cell(n_net, n_ratio);
for ni=1:n_net
    for ri=1:n_ratio
        res_name = sprintf('./SyNet_OverlapPerThresh/%sOV-PerThresh_%s_%s_MP%d_%s_RV%0.2f.mat', Overlap_Check, Ref_Name, Net_lst{ni}, MAX_N_PAIR, Limit_method, Ratio_lst(ri));
        %res_name = sprintf('./SyNet_OverlapPerThresh/LnkOV-PerThresh_%s_%s_MP%d_%s_RV%0.2f.mat', Ref_Name, Net_lst{ni}, MAX_N_PAIR, Limit_method, Ratio_lst(ri));
        fprintf('Loading [%s]\n', res_name);
        
        if strcmp(Overlap_Check, 'Gene')
            res_data = load(res_name, 'Real_Freq', 'Rand_Freq', 'n_SelGene');
            param_str{ni, ri} = sprintf('%d/%d', res_data.Real_Freq, res_data.n_SelGene);
        else
            res_data = load(res_name, 'Real_Freq', 'Rand_Freq', 'n_TopPair');
            param_str{ni, ri} = sprintf('%d/%d', res_data.Real_Freq, res_data.n_TopPair);
        end
        
        if std(res_data.Rand_Freq)<0.001
            continue;
        end
        
        zscr_Mat(ni, ri) = (res_data.Real_Freq - mean(res_data.Rand_Freq)) / std(res_data.Rand_Freq);
    end
end


%% Plotting
figure('Position', [100 100 1400 500]);
hold on
for ni=1:n_net
    [clr_map, Method_lst{ni,1}] = getColor(Net_lst{ni});
    plot(1:n_ratio, zscr_Mat(ni,:), [Mrk_lst{ni} '-'], 'Color', clr_map, 'MarkerSize', 10, 'MarkerFaceColor', clr_map);
    
    [~, best_ind] = max(zscr_Mat(ni,:));
    text(best_ind, zscr_Mat(ni, best_ind), param_str{ni, best_ind}, ...
        'VerticalAlignment', 'Bottom', 'HorizontalAlignment', 'Center', 'FontWeight', 'Bold');
end
set(gca, 'XTick', 1:n_ratio, 'XTickLabel', Ratio_lst*100);
xlim([0 n_ratio+1]);
xlabel('Ratio of links');
ylabel('Z-score');
title(sprintf('Overlap of [%s] and biological networks [Over %s, Method = %s]', Ref_Name, Overlap_Check, Limit_method));
legend(Method_lst);


