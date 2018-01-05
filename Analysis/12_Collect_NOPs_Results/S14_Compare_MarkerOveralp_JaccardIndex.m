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
y_lim = [0.6 0.66];

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
    [met_clr, Method_lbl{si,1}] = getColor(net_info{1});
    
    index_freq = accumarray(nonzeros(res_data.Marker_lst(:)), 1);
    [mrk_frq, mrk_sind] = sort(index_freq, 'Descend');
    n_mrk = round(median(sum(res_data.Marker_lst~=0,2)));
    Top_Mrk{si} = mrk_sind(1:n_mrk);
end

%% Check overlap between markers
Jac_Mat = nan(n_res);
for si=1:n_res
    for sj=1:si-1
        uni_set = union(Top_Mrk{si}, Top_Mrk{sj});
        int_set = intersect(Top_Mrk{si}, Top_Mrk{sj});
        Jac_Mat(si, sj) = numel(int_set) / numel(uni_set);
    end
end

%% Plot jaccard matrix
img_h = imagesc(Jac_Mat);
set(img_h, 'AlphaData', ~isnan(Jac_Mat));
set(gca, 'FontSize', 12, 'FontWeight', 'Bold', 'Clim', [0 0.4], 'XTick', 1:n_res, 'XTickLabel', Method_lbl, ...
    'YTick', 1:n_res, 'YTickLabel', Method_lbl, 'Box', 'off');
colormap(jet(10));
clrbar_h = colorbar('North');
xlabel(clrbar_h, 'Jaccard Index');
view(-45, 90);

%% Saving
output_name = sprintf('./Plots/S14_CompareOverlap_BetweenIdentifiedMarkers_GenesInNet.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [9 7], 'PaperPosition', [0 0 9 7]);
print('-dpdf', '-r300', output_name);




