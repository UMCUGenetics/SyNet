clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/Tight_Subplot/');
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
neg_clr = [0.7 0.7 0.7];
net_lst = {'HumanInt' 'BioPlex','BioGRID','IntAct','STRING','HBBrain','HBKidney','HBOvary','HBLympNode','HBGland'};
n_net = numel(net_lst);
tm_lst = {'Degree' 'PageRank-FB0.85' 'Betweenness'}; 
n_tm = numel(tm_lst);

%% Loading TM data
load('./Topological_Data/TMData_NS20000_NF180.mat');
[n_sample, n_feature] = size(TM_Data_filtered);

%% Main loop
figure('Position', [100 100 1600 500]);
sp_h = tight_subplot(1, 3, 0.03, 0.1, 0.02);
step = 1;
for ti=1:n_tm
    set(gcf, 'CurrentAxes', sp_h(step));
    hold on
    Method_lst = {};
    for ni=1:n_net
        Feat_Ind = cellfun('length', strfind(TM_Name, sprintf('%s-%s-A', net_lst{ni}, tm_lst{ti})))==1;
        if sum(Feat_Ind)~=1, error(); end
        
        data_neg = TM_Data_filtered(TM_Label==-1, Feat_Ind);
        box_h = BoxPlotEx(data_neg, 'Positions', ni-0.17, 'Color', neg_clr, 'Symbol', '', 'Widths', 0.3);
        set(box_h, 'LineWidth', 2);
        
        data_pos = TM_Data_filtered(TM_Label==1, Feat_Ind);
        [met_clr, Method_lst{ni,1}] = getColor(net_lst{ni});
        box_h = BoxPlotEx(data_pos, 'Positions', ni+0.17, 'Color', met_clr, 'Symbol', '', 'Widths', 0.3);
        set(box_h, 'LineWidth', 2);        
    end
    switch tm_lst{ti}
        case 'Degree'
            set(gca, 'YLim', [0 250]);
        case 'PageRank-FB0.85'
            set(gca, 'YLim', [0 1.5e-3]);
        case 'Eigenvector'
            set(gca, 'YLim', [0 1e-3]);
        case 'Closeness'
            set(gca, 'YLim', [0 2.3e-5]);
        case 'Betweenness'
            set(gca, 'YLim', [0 4e4]);
    end
    set(gca, 'XTick', 1:n_net, 'XTickLabel', Method_lst, 'XTickLabelRotation', 20, ...
        'XLim', [0.3 n_net+0.7], 'FontWeight', 'Bold');
    title(tm_lst{ti}, 'FontSize', 14);
    step = step + 1;
end

% return
%% Saving
output_name = sprintf('./Plots/S03_Enrichment_Evaluation_SyNet_Over_TM.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [15 4], 'PaperPosition', [0 0 15 4]);
print('-dpdf', '-r300', output_name);

