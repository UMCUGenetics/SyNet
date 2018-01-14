clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/getTop/');

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
n_study = 14;
n_rep = 10;

%% Collecting top markers
Jac_Dist = cell(n_net, 1);
for ni=1:n_net
    res_name = [res_path net_lst{ni}];
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    if size(res_data.Marker_lst,1)~=n_epoch, error(); end
    
    res_info = regexp(net_lst{ni}, '_', 'split');
    net_info = regexp(res_info{4}, '-', 'Split');
    [met_clr(ni,:), Method_lbl{ni,1}] = getColor(net_info{1});
    
    ei = 1;
    mrk_perEpoch = cell(n_study, n_rep);
    mrk_set = {};
    for si=1:n_study
        for ri=1:n_rep
            %mrk_perEpoch{si, ri} = nonzeros(res_data.Marker_lst(ei, :));
            mrk_set{end+1,1} = nonzeros(res_data.Marker_lst(ei, :));
            ei = ei + 1;
        end
        %mrk_set{end+1,1} = vertcat(mrk_perEpoch{si,:});
        %[mrk_lst, mrk_freq] = getTop(vertcat(mrk_perEpoch{si,:}));
        %mrk_set{end+1,1} = mrk_lst(mrk_freq>2);
        %mrk_set{end+1,1} = getTop(vertcat(mrk_perEpoch{si,:}), 100);
    end
    
    n_batch = numel(mrk_set);
    jac_index = cell(n_batch, 1);
    for bi=1:n_batch
        for bj=bi+1:n_batch
            uni_set = union(mrk_set{bi}, mrk_set{bj});
            int_set = intersect(mrk_set{bi}, mrk_set{bj});
            jac_index{bi} = [jac_index{bi}; numel(int_set) / numel(uni_set)];
        end
    end
    Jac_Dist{ni} = vertcat(jac_index{:});
end

%% Plot jaccard matrix
close all
figure('Position', [100 100 800 600]);
hold on
% figure('Position', [100 100 1500 400]);
% hold on
% for ni=1:n_net
%     BoxPlotEx(Jac_Dist{ni}, 'Color', met_clr(ni,:), 'Positions', ni, 'Widths', 0.3, 'Symbol', ''); %
%     text(ni, 0.2, sprintf('%0.3f', std(Jac_Dist{ni})));
% end
char_lst = {'O' 'x' '>' '<' '^' 'v' '+' 'o' 'd' 's' '^' 'p' 'x'};
for ni=1:n_net
    if ismember(char_lst{ni}, {'x' '+'})
        plt_h(ni,1) = plot(mean(Jac_Dist{ni}), std(Jac_Dist{ni}), char_lst{ni}, ...
            'MarkerSize', 10, ...
            'MarkerFaceColor', met_clr(ni,:), 'MarkerEdgeColor', met_clr(ni,:), 'LineWidth', 3);
    elseif ismember(char_lst{ni}, {'p'})
        plt_h(ni,1) = plot(mean(Jac_Dist{ni}), std(Jac_Dist{ni}), char_lst{ni}, ...
            'MarkerSize', 13, ...
            'MarkerFaceColor', met_clr(ni,:), 'MarkerEdgeColor', 'none');
    else
        plt_h(ni,1) = plot(mean(Jac_Dist{ni}), std(Jac_Dist{ni}), char_lst{ni}, ...
            'MarkerSize', 8, ...
            'MarkerFaceColor', met_clr(ni,:), 'MarkerEdgeColor', 'none');
    end
end
legend(Method_lbl);
xlabel('Average overlap (Jaccard index)');
ylabel('Standard deviation of overlap (Jaccard index)');
set(gca, 'FontWeight', 'Bold', 'FontSize', 10);
title('Overlap between identified markers', 'FontSize', 12, 'FontWeight', 'Bold');

% return
% %% Saving
output_name = sprintf('./Plots/S16_CompareOverlap_BetweenIdentifiedMarkers_AcrossFolds_ScatterPlot.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [6 5], 'PaperPosition', [0 0 6 5]);
print('-dpdf', '-r300', output_name);




