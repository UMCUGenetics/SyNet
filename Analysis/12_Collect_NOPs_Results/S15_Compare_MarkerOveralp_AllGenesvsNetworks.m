clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    'MRK_CVT01_NetLasso_HumanInt-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_BioPlex-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_BioGRID-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_IntAct-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_STRING-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBOvary-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBBrain-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBKidney-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBGland-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_HBLympNode-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_ACr-P50000_MSN-500.mat'
    'MRK_CVT01_NetLasso_AvgSynACr-P50000_MSN-500.mat'
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
    };
n_res = numel(res_lst);

%% Collecting top markers
Top_Mrk = cell(n_res, 1);
for si=1:n_res
    res_name = [res_path res_lst{si}];
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    auc_rep = mean(res_data.AUC_mat, 1);
    res_info = regexp(res_lst{si}, '_', 'split');
    net_info = regexp(res_info{4}, '-', 'Split');
    [met_clr(si,:), Method_lbl{si,1}] = getColor(net_info{1});
    
    index_freq = accumarray(nonzeros(res_data.Marker_lst(:)), 1);
    [mrk_frq, mrk_sind] = sort(index_freq, 'Descend');
    n_mrk = round(median(sum(res_data.Marker_lst~=0,2)));
    Top_Mrk{si} = mrk_sind(1:n_mrk);
end

%% Plot jaccard matrix
close all
figure('Position', [100 100 1500 400]);
hold on
for si=1:n_res-1
    uni_set = union(Top_Mrk{end}, Top_Mrk{si});
    int_set = intersect(Top_Mrk{end}, Top_Mrk{si});
    Jac_Index = numel(int_set) / numel(uni_set);
    
    bar(si, Jac_Index, 'FaceColor', met_clr(si,:), 'BarWidth', 0.7, 'EdgeColor', met_clr(si,:)*0.6);
end

xlim([0 n_res]);
ylim([0 0.5]);
set(gca, 'XTick', 1:n_res-1, 'XTickLabel', Method_lbl(1:n_res-1), ...
    'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('Jaccard Index', 'FontWeight', 'Bold');
title('Overlap comparison between identified markers (all genes vs. networks)', 'FontSize', 12, 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S14_CompareOverlap_BetweenIdentifiedMarkers_AllvsInNet.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [14 4], 'PaperPosition', [0 0 14 4]);
print('-dpdf', '-r300', output_name);




