clc;
clear;

%% 
addpath('../_Utilities/');
TP_Data = load('./Top_Pairs/Top_SyNet.mat');
Pair_Info = TP_Data.PP_Info;

%% Calculate axes
% Pair_Info(:,6) = mean(Pair_Info(:,3:4), 2);
Pair_Info(:,6) = min(Pair_Info(:,3:4), [], 2);

%% Compute score
ox = Pair_Info(:,6);
ox = ox-min(ox);
ox = ox/max(ox);

oy = Pair_Info(:,7);
oy = oy-min(oy);
oy = oy/max(oy);

GE_data = load(getPath('SyNet'), 'Gene_Expression');
cr_mat = abs(corr(GE_data.Gene_Expression, 'Type', 'Spearman'));
n_gene = size(cr_mat,1);
pind = sub2ind([n_gene n_gene], Pair_Info(:,1), Pair_Info(:,2));
Pair_Info(:,9) = cr_mat(pind);
oz = cr_mat(pind);
oz = oz-min(oz);
oz = oz/max(oz);

Pair_Info(:,10) = -sqrt((1-ox).^2 + (1-oy).^2 + (1-oz).^2);

%% Sorting
[~, sid] = sort(Pair_Info(:,10), 'Descend');
PP_Info = Pair_Info(sid, :);

%% Top pairs
TP_Data.Gene_Name(PP_Info(1:100,1:2))

%% Gather Top genes
Top_GInd = unique(PP_Info(:,1:2)', 'Stable');
GN = TP_Data.Gene_Name(Top_GInd(1:500));

%% Visualize
close all
plot(PP_Info(:,6), PP_Info(:,7), '.', 'Color', [0.8 0.8 0.8]);
hold on
PP_Info(10000:end, :) = [];
plot(PP_Info(:,6), PP_Info(:,7), '.', 'Color', [0.2 0.9 0.2]);


%% Choose top genes
Gene_Info = [[Pair_Info(:,1); Pair_Info(:,2)] [Pair_Info(:,3); Pair_Info(:,4)]];
[g_val, g_ind] = sort(Gene_Info(:,2), 'Descend');
Gene_Info = unique(Gene_Info(g_ind,:), 'rows', 'Stable');
GN = TP_Data.Gene_Name(Gene_Info(1:500,1));