clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
ge_name = 'SyNet';

%% Load top pairs
load(['./Top_Pairs/Top_' ge_name '.mat'], 'PP_Info', 'Gene_Name');
n_pairs = size(PP_Info, 1);
n_top = 500;

%% Take top pairs based on observation
Top_Gene = unique(PP_Info(:,1:2)', 'Stable');
Top_Gene(n_top+1:end) = [];
Top_Gene = flipud(Top_Gene);
Pair_Score = zeros(n_pairs, 1);
for ti=1:n_top
	is_in = PP_Info(:,1)==Top_Gene(ti) | PP_Info(:,2)==Top_Gene(ti);
	Pair_Score(is_in) = ti;
end

%% Plotting
close all
figure();
clr_map = [0.8 0.8 0.8; jet(n_top)];
%clr_map = AdvancedColormap('g');
set(gca, 'FontWeight', 'Bold');
xlabel('Pairwise AUC');
ylabel('Synergy');
title(['Score Space (' ge_name  ')']);
hold on
for ci=0:n_top
	is_cls = find(Pair_Score==ci);
	plot(PP_Info(is_cls,6), PP_Info(is_cls,7), '.b', 'MarkerEdgeColor', clr_map(ci+1,:), 'MarkerSize', 1);
	is_cls(is_cls>10000) = [];
	plot(PP_Info(is_cls,6), PP_Info(is_cls,7), '+b', 'MarkerEdgeColor', clr_map(ci+1,:), 'MarkerSize', 2);
end
colormap(clr_map);
clr_h = colorbar();
caxis([0 n_top]);
set(clr_h, 'YTick', 0:50:n_top);
% xlim([0.58 0.67]);
% ylim([0.99 1.1]);
ylabel(clr_h, 'Score');

%% Saving the plot
sav_name = ['./Plots/ScoreSpace_' ge_name '_BasedOnObserved.png'];
print(gcf, '-dpng', '-r300', sav_name);

