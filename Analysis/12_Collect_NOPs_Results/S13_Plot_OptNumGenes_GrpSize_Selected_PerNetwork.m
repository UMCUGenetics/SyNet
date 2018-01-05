clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/Tight_Subplot/');

%% Select methods and networks
res_path = './Collected_Results/';
res_lst = {
    'MRK_CVT01_NetGL_HumanInt-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_BioPlex-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_BioGRID-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_IntAct-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_STRING-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_HBOvary-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_HBBrain-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_HBKidney-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_HBGland-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_HBLympNode-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_ACr-P50000_MSN-500.mat'
    'MRK_CVT01_NetGL_AvgSynACr-P50000_MSN-500.mat'
    };
n_res = numel(res_lst);
Param_lst = {
    [50 100 300 500 700 1500 3000]
    [2 3 5 7 10]
    };
Param_FName = {'OptNGene' 'OptGrpSize'};
Param_yLabel = {'Number of genes' 'Group size'};
n_param = numel(Param_lst);

%% Collect data
figure('Position', [100 100 1400 600]);
sp_h = tight_subplot(1, 2, 0.05, 0.1, 0.1, 'TickLength', [0.0 0.0]);
for oi=1:n_param
    n_item = numel(Param_lst{oi});
    Opt_mat = zeros(n_item, n_res);
    Net_Name = cell(n_res, 1);
    for ni=1:n_res
        res_name = sprintf([res_path '%s'], res_lst{ni});
        fprintf('Reading [%s]\n', res_name);
        net_data = load(res_name);
        if any(isnan(net_data.AUC_mat(:))), error(); end
        res_info = regexp(res_lst{ni}, '[_-]', 'Split');
        [~, Net_Name{ni}] = getColor(res_info{4});
        
        for ii=1:n_item
            is_in = net_data.(Param_FName{oi})==Param_lst{oi}(ii);
            Opt_mat(ii, ni) = sum(is_in);
        end
        if sum(Opt_mat(:,ni))~=140, error(); end
    end
    
    %% Plot matrix
    set(gcf, 'CurrentAxes', sp_h(oi));
    imagesc(Opt_mat);
    colormap(parula(8));
    clr_h = colorbar('XTick', 0:5:60);
    ylabel(clr_h, 'Frequency', 'FontWeight', 'Bold');
    set(gca, 'XTick', 1:n_res, 'XTickLabel', Net_Name, 'FontWeight', 'Bold', 'TickLength', [0.0 0.0], ...
        'YTick', 1:n_item, 'YTickLabel', Param_lst{oi}, 'CLim', [0 50]);
    ylabel(Param_yLabel{oi});
    title(sprintf('Optimal %s', lower(Param_yLabel{oi})));
end

% return
%% Saving
output_name = sprintf('./Plots/S13_OptimalNumber_Genes-GroupSize.pdf');
set(gcf, 'PaperUnit', 'inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [22 4], 'PaperPosition', [0 0 22 4]);
print('-dpdf', '-r300', output_name);

