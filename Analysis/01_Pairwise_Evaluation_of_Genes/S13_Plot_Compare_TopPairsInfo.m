clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
ge_name = 'SyNet';

%% Load DSN
% net_name = ['./Network_Files/DSN_' ge_name '.mat'];
net_name = ['./Network_Files/DSN_' ge_name 'S01.mat'];
fprintf('Loading DSN [%s]\n', net_name);
dsn_data = load(net_name, 'Pair_AUC', 'Pair_Std', 'Gene_Name');
n_gene = size(dsn_data.Pair_AUC, 1);
n_total = n_gene*(n_gene-1)/2;
fprintf('In total [%d] gene pairs exist.\n', n_total);

%% Identify pairs
fprintf('Loading pair list.\n');
Pair_Info = zeros(n_total, 10, 'single');
[Pair_Info(:,1), Pair_Info(:,2)] = find(triu(ones(n_gene), 1));

%% Load top pairs
% tp_data = load(['./Top_Pairs/Top_' ge_name '.mat'], 'PP_Info', 'Gene_Name');

%% Add scores
%% Generate Pair Matrix
Ind_AUC = dsn_data.Pair_AUC(1:n_gene+1:end)';
Pair_Info(:,3) = Ind_AUC(Pair_Info(:,1));
Pair_Info(:,4) = Ind_AUC(Pair_Info(:,2));
Pair_Info(:,5) = max(Pair_Info(:,3:4), [], 2);
pair_ind = sub2ind([n_gene n_gene], Pair_Info(:,1), Pair_Info(:,2));
Pair_Info(:,6) = dsn_data.Pair_AUC(pair_ind);
Pair_Info(:,7) = Pair_Info(:,6)./Pair_Info(:,5);
%clear Pair_AUC

%% Normalizing scores
ox = Pair_Info(:,6);
ox = ox-min(ox);
ox = ox/max(ox);
oy = Pair_Info(:,7);
oy = oy-min(oy);
oy = oy/max(oy);
Pair_Info(:,8) = -sqrt((1-ox).^2 + (1-oy).^2);
Pair_Info(:,9) = dsn_data.Pair_Std(pair_ind);

%% Add correlation
ge_path = getPath(ge_name);
ge_data = load(ge_path, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
crr_mat = corr(ge_data.Gene_Expression, 'Type', 'Spearman');
Pair_Info(:,10) = crr_mat(pair_ind);

%% Get Top pairs
[~, sid] = sort(-Pair_Info(:,9), 'Descend');
% Top_Pair = Pair_Info(sid(1:10000), :);
Top_Pair = Pair_Info(sid, :);

%% Discretize
Top_Pair(Top_Pair(:,7)<1.05,:) = [];
Top_Pair(abs(Top_Pair(:,10))<0.05,:) = [];
clf
n_xbin = 10;
n_ybin = 5;
[Freq_mat, Xedges, Yedges, XClass, YClass] = plotAxisCvg(Top_Pair(:, [6 8]), n_xbin, n_ybin);
img = Freq_mat./sum(Freq_mat,2);
imagesc(img');
set(gca, 'XTick', 1:n_xbin, 'XTickLabel', round(Xedges,3), 'XTickLabelRotation', 15, ...
	'YTick', 1:n_ybin, 'YTickLabel', round(Yedges,3), 'YDir', 'Normal');

