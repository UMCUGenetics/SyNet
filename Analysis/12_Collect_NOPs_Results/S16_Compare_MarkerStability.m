clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');

%% Select methods and networks
res_path = './Collected_Results/';
net_lst = {
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
n_net = numel(net_lst);
n_epoch = 140;
n_TopMrk = 25;
n_study = 14;
n_rep = 10;

%% Collecting top markers
Jac_Dist = cell(n_net, 1);
for ni=1:n_net
    res_name = [res_path net_lst{ni}];
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    if size(res_data.Marker_lst,1)~=140, error(); end
    
    res_info = regexp(net_lst{ni}, '_', 'split');
    net_info = regexp(res_info{4}, '-', 'Split');
    [met_clr(ni,:), Method_lbl{ni,1}] = getColor(net_info{1});
    
    ei = 1;
    mrk_set = cell(n_study, n_rep);
    for si=1:n_study
        for ri=1:n_rep
            mrk_set{si,ri} = res_data.Marker_lst(ei, 1:n_TopMrk);
            ei = ei + 1;
        end
    end
    
    jac_index = cell(n_study, n_rep);
    jac_med = zeros(n_study, n_rep);
    for si=1:n_study
        for ri=1:n_rep
            for sj=1:n_study
                for rj=1:n_rep
                    if si==sj && ri==rj, continue; end
                    uni_set = union(mrk_set{si,ri}, mrk_set{sj,rj});
                    int_set = intersect(mrk_set{si,ri}, mrk_set{sj,rj});
                    jac_index{si,ri} = [jac_index{si,ri}; numel(int_set) / numel(uni_set)];
                end
            end
            jac_med(si, ri) = median(jac_index{si,ri});
        end
    end
    
    Jac_Dist{ni} = jac_med(:);
    %Jac_Dist{ni} = vertcat(jac_index{:});
end

%% Plot jaccard matrix
close all
figure('Position', [100 100 1500 400]);
hold on
for ni=1:n_net
    BoxPlotEx(Jac_Dist{ni}, 'Color', met_clr(ni,:), 'Positions', ni, 'Widths', 0.3); % , 'Symbol', ''
end

xlim([0 n_net+1]);
ylim([0 0.3]);
set(gca, 'XTick', 1:n_net, 'XTickLabel', Method_lbl(1:n_net), ...
    'FontWeight', 'Bold', 'FontSize', 10, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('Jaccard Index', 'FontWeight', 'Bold');
title('Overlap between identified markers across folds', 'FontSize', 12, 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S16_CompareOverlap_BetweenIdentifiedMarkers_AcrossFolds.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [14 4], 'PaperPosition', [0 0 14 4]);
print('-dpdf', '-r300', output_name);




