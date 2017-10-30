clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
res_name = './Collected_Results/MRK_CVT01_LExAG_None-G11748_MSN-500.mat';

%% Loading the results
fprintf('Loading results from [%s] ...\n', res_name);
clc_data = load(res_name);
AUC_mat = clc_data.AUC_mat;
n_study = numel(Study_Name);

%% Sort studies
[~, sid] = sort(mean(AUC_mat,2));
AUC_mat = AUC_mat(sid,:);
Study_Name = Study_Name(sid);

%% Box plot
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(n_study);
X_Label = cell(n_study,1);
for si=1:n_study
    val = AUC_mat(si,:)';
    if any(isnan(val)), error(); end
    Point_Param.Colormap = clr_map(si,:);
    Point_Param.ColorCAxis = [min(val) max(val)];
    box_h = boxplotEx(val, si, {}, Point_Param);
    set(box_h, 'Color', clr_map(si,:), 'Marker', 'none');
    X_Label{si} = sprintf('%s', Study_Name{si});
end
X_Label = regexprep(X_Label, 'ACES;', '');
set(gca, 'XLim', [0 n_study+1], 'YLim', [0.50 0.81]);
set(gca, 'XTick', 1:n_study, 'XTickLabel', X_Label, 'XTickLabelRotation', 15, 'FontWeight', 'Bold', 'FontSize', 12);
ylabel('AUC', 'FontSize', 12, 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S10_AUC-Performance_CrossStudy.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [15 4], 'PaperPosition', [0 0 15 4]);
print('-dpdf', '-r300', output_name);

