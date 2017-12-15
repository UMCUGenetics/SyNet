clc;
clear;
close all

%% Initialization
% addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Distribution_Plot/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/Hatchfill/');
addpath('../_Utilities/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Select methods and networks
res_lst = {
%     {'MRK_CVT01_NetLasso_BioPlex-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_BioPlex-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_BioGRID-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_BioGRID-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_IntAct-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_IntAct-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_STRING-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_STRING-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_HBOvary-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBOvary-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_HBBone-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBBone-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBGland-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBGland-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_HBLympNode-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBLympNode-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_HBEpith-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBEpith-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_ACr-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_ACr-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_HBColon-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBColon-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_HBBrain-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBBrain-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_HBIntestine-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBIntestine-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_HBLung-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_HBLung-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_Syn-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_Syn-P25000_MSN-500.mat'}
%     {'MRK_CVT01_NetLasso_AvgSyn-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_AvgSyn-P25000_MSN-500.mat'}
    {'MRK_CVT01_NetLasso_AvgSynACr-P25000_MSN-500.mat' 'MRK_CVT01_NetGL_AvgSynACr-P25000_MSN-500.mat'}
    % {'MRK_CVT01_Lasso_AvgSynACr-P10000_MSN-500.mat' 'MRK_CVT01_GLasso_AvgSynACr-P10000_MSN-500.mat'}
    {'MRK_CVT01_LExAG_None-G11748_MSN-500.mat' 'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'}
    };
n_res = numel(res_lst);
y_lim = [0.61 0.67];

%% Plotting performance
close all
figure('Position', [100 100 1500 400]);
hold on
X_lbl = {};
xlim([0.5 n_res+0.5]);
ylim(y_lim);
for si=1:n_res
    auc1_study = LoadData(res_lst{si}{1});
    auc2_study = LoadData(res_lst{si}{2});
    res_info = regexp(res_lst{si}{2}, '_', 'split');
    net_info = regexp(res_info{4}, '-', 'Split');
    X_lbl{si,1} = sprintf('%s', net_info{1});
    [met_clr, met_lbl] = getColor(X_lbl{si});
    
    bar_pos = si-0.2;
    bar_h = patch(bar_pos+[-1 1 1 -1]*0.18, [0 0 mean(auc1_study) mean(auc1_study)], 'b', 'FaceColor', met_clr);
    errorbarEx(bar_pos, mean(auc1_study), std(auc1_study), std(auc1_study), 1.5, 0.1, [0 0 0]);
    
    if si~=n_res
        bar_pos = si+0.2;
        bar_h = patch(bar_pos+[-1 1 1 -1]*0.18, [0 0 mean(auc2_study) mean(auc2_study)], 'b', 'FaceColor', met_clr);
        patch_h = findobj(bar_h, 'Type', 'patch');
        hatch_h = hatchfill(patch_h, 'single', -45, 10, min(met_clr*1.1, [1 1 1]));
        %set(hatch_h, 'Color', met_clr*0.8, 'LineWidth', 1.5);
        errorbarEx(bar_pos, mean(auc2_study), std(auc2_study), std(auc2_study), 1.5, 0.1, [0 0 0]);
    end
    text(si, y_lim(1)-0.002, met_lbl, 'FontSize', 10, 'FontWeight', 'Bold', 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top', 'Rotation', 0);
end

y_tick = get(gca, 'YTick');
Y_lbl = arrayfun(@(y) sprintf('%0.0f%%', y*100), y_tick, 'UniformOutput', 0);
set(gca, 'XTick', [], 'XTickLabel', [], 'XTickLabelRotation', 10, ...
    'YTick', y_tick, 'YTickLabel', Y_lbl, 'FontWeight', 'Bold', 'FontSize', 7, ...
    'Ygrid', 'on', 'GridColor', [0.7 0.7 0.7], 'GridAlpha', 0.4);
ylabel('AUC', 'FontWeight', 'Bold');

return
%% Saving
output_name = sprintf('./Plots/S06_PerformanceComparison_NetBased.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [5 3], 'PaperPosition', [0 0 5 3]);
print('-dpdf', '-r300', output_name);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function auc_study = LoadData(res_name)
res_path = './Collected_Results/';

res_name = [res_path res_name];
fprintf('Reading [%s]\n', res_name);
res_data = load(res_name);
if any(isnan(res_data.AUC_mat(:))), error(); end
auc_study = mean(res_data.AUC_mat,1);
end
