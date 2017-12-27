clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');

%% Select methods and networks
res_path = './Collected_Results/';
Net_lst = {
    'HumanInt'
    'BioPlex'
    'BioGRID'
    'IntAct'
    'STRING'
    'HBBrain'
    'HBKidney'
    'HBOvary'
    'HBGland'
    'HBLympNode'
    'ACr'
    'AvgSynACr'
    };
n_net = numel(Net_lst);
Param_lst = {
    [50 100 300 500 700 1500 3000]
    [2 3 5 7 10]
    };
Param_FName = {'OptNGene' 'OptGrpSize'};
n_param = numel(Param_lst);

%% Collect data
figure('Position', [100 100 1400 600]);
for oi=1:n_param
    n_item = numel(Param_lst{oi});
    Opt_mat = zeros(n_item, n_net);
    Net_Name = cell(n_net, 1);
    for ni=1:n_net
        res_name = sprintf([res_path 'MRK_CVT01_NetGL_%s-P50000_MSN-500.mat'], Net_lst{ni});
        fprintf('Reading [%s]\n', res_name);
        net_data = load(res_name);
        if any(isnan(net_data.AUC_mat(:))), error(); end
        [~, Net_Name{ni}] = getColor(Net_lst{ni});
        
        for ii=1:n_item
            is_in = net_data.(Param_FName{oi})==Param_lst{oi}(ii);
            Opt_mat(ii, ni) = sum(is_in);
        end
        if sum(Opt_mat(:,ni))~=140, error(); end
    end
    
    %% Plot matrix
    subplot(1, n_param, oi);
    imagesc(Opt_mat);
    % colormap(bone(5));
    set(gca, 'XTick', 1:n_net, 'XTickLabel', Net_Name, ...
        'YTick', 1:n_item, 'YTickLabel', Param_lst{oi});
    title(Param_FName{oi});
end

return
%% Saving
output_name = sprintf('./Plots/S05_04_NumOptGenes.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [13 3], 'PaperPosition', [0 0 13 3]);
print('-dpdf', '-r300', output_name);

