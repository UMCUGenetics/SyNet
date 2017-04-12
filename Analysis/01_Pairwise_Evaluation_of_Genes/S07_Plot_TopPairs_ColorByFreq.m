clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
ge_name = 'SyNet';

%% Load top pairs
load(['./Top_Pairs/Top_' ge_name '.mat'], 'PP_Info', 'Gene_Name');
n_top = 10000;

%% Set labels
TG_Name = {};
for ti=1:n_top
	if ~ismember(Gene_Name{PP_Info(ti,1)}, TG_Name)
		TG_Name{end+1,1} = Gene_Name{PP_Info(ti,1)};
	end
	if numel(TG_Name)>=100
		break;
	end
	if ~ismember(Gene_Name{PP_Info(ti,2)}, TG_Name)
		TG_Name{end+1,1} = Gene_Name{PP_Info(ti,2)};
	end
	if numel(TG_Name)>=100
		break;
	end
end

Pair_Rank = zeros(n_top, 1);
[Freq_item, Freq_freq] = getTop(PP_Info(1:n_top,1:2), 100);
Freq_item = flipud(Freq_item);
Freq_freq = flipud(Freq_freq);
n_freq = numel(Freq_item);
for fi=1:n_freq
	is_in = PP_Info(:,1)==Freq_item(fi) | PP_Info(:,2)==Freq_item(fi);
	Pair_Rank(is_in) = fi;
end

%% Plotting
close all
figure();
clr_map = [0.8 0.8 0.8; jet(n_freq)];
%clr_map = AdvancedColormap('g');
set(gca, 'FontWeight', 'Bold');
xlabel('Pairwise AUC');
ylabel('Synergy');
title(['Score Space (' ge_name  ')']);
hold on
for ci=0:n_freq
	is_cls = find(Pair_Rank==ci);
	plot(PP_Info(is_cls,6), PP_Info(is_cls,7), '.b', 'MarkerEdgeColor', clr_map(ci+1,:), 'MarkerSize', 1);
	is_cls(is_cls>10000) = [];
	plot(PP_Info(is_cls,6), PP_Info(is_cls,7), '+b', 'MarkerEdgeColor', clr_map(ci+1,:), 'MarkerSize', 2);
end
colormap(clr_map);
clr_h = colorbar();
caxis([0 n_freq]);
set(clr_h, 'YTick', [1:5:n_freq 100], 'YTickLabel', [Freq_freq(1:5:end); Freq_freq(end)]);
% xlim([0.58 0.67]);
% ylim([0.99 1.1]);
ylabel(clr_h, 'Frequency');

%% Saving the plot
sav_name = ['./Plots/ScoreSpace_' ge_name '.png'];
print(gcf, '-dpng', '-r300', sav_name);
