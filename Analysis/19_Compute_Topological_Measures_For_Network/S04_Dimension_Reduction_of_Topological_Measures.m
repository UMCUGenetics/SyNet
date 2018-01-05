clc;
clear;

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
Perp_lst = [5 10 20 40 70 110 200 400];
mrk_lst = {'+' 'o'};
n_mrk = numel(mrk_lst);

%% Load TM data
load('./Topological_Data/TMData-OneGRND_NS7088_NF180.mat', 'TM_Data_z', 'TM_Label', 'TM_Name');
qTM_Data = quantilenorm(TM_Data_z); % , 'Display', true

%% Reducing dimension Using PCA
% [coeff,score,latent,tsquared,explained,mu] = pca(zTM_Data, 'NumComponents', 40);

%% Reducing dimension Using t-SNE
for pi=1:numel(Perp_lst)
    perplexity = Perp_lst(pi);
    fprintf('Perplexity set to: %d\n', perplexity);
    initial_dims = 30;
    no_dims = 2;
    pData = fast_tsne(qTM_Data, no_dims, initial_dims, perplexity, 0.4, 500);
    pLabel = (TM_Label==1)+1;
    % theta = 0 corresponds to standard, slow t-SNE, while theta = 1 makes very crude approximations.
    
    %% Plot the reduced data
    close all
    clr_map = [0.7 0.7 0.7; 0.3 0.8 0.2];
%     figure('Visible', 'off');
    hold on
    plot_h = [];
    for ci=1:2
        plot_h(ci,1) = scatter(pData(pLabel==ci,1), pData(pLabel==ci,2), 10, mrk_lst{ci}, 'MarkerFaceColor', clr_map(ci,:), 'MarkerEdgeColor', clr_map(ci,:), 'MarkerFaceAlpha', 0.7);
    end
    set(plot_h(2), 'MarkerEdgeColor', 'None');
    legend(plot_h([2 1]), {'SyNet' 'Random'}, 'fontsize', 12, 'fontweight', 'bold');
    title('Visualization of topological measures using t-SNE');
    
    %% Saving
    sav_name = sprintf('./Plots/S04_DR_TopologicalMeasure_t-SNE_Perp%03d.png', perplexity);
    print(gcf, '-dpng', '-r300', sav_name);
end


