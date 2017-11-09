clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
S04_Compare_AUC_Stability;
ge_path = getPath('ACES');
MAX_N_SUBNET = 50;
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
Target_lst = strrep(Study_Name, 'ACES;', '');

%% Load data
fprintf('Loading data: [%s]\n', ge_path);
load(ge_path, 'Gene_Expression', 'Patient_Label', 'Gene_Name', 'Study_Name', 'Study_Index');
n_study = max(Study_Index);
zData = zscore(Gene_Expression);
n_gene = size(zData, 2);
n_feat = min([MAX_N_SUBNET n_gene]);

%% Main loop over studies
te_auc = zeros(n_study, 1);
for si=1:n_study
	tar_ind = find(ismember(Study_Name, Target_lst{si}));
	iTr = Study_Index ~= tar_ind;
	iTe = Study_Index == tar_ind;
	
	lTr = (Patient_Label(iTr)==1)*2-1;
	lTe = (Patient_Label(iTe)==1)*2-1;
	
	zTr = zData(iTr,:);
	zTe = zData(iTe,:);
	
	%% Select top genes
	fprintf('Evaluting [%d] individual genes.\n', n_gene);
	pv_vec = ttest2Ex(zTr, lTr);
    
	%% Selecting top genes
	[~, scr_ind] = sort(-log10(pv_vec), 'Descend');
	sTr = zTr(:, scr_ind(1:n_feat));
	sTe = zTe(:, scr_ind(1:n_feat));
	
	%% Traning the final model
	fprintf('Training the final model over [%d] features...\n', n_feat);
	pred = nmc(sTr, lTr, sTe);
	
	%% Saving results
	te_auc(si,1) = getAUC(lTe, pred, 50);
end

%% Plot the results
bar_h = bar(2:2:n_study*2, te_auc, 'g', 'BarWidth', 0.4);
bar_h.FaceAlpha = 0.4;
bar_h.EdgeColor = 'none';
set(gca, 'XTick', 2:2:n_study*2, 'XTickLabel', Target_lst, 'XTickLabelRotation', 20);

