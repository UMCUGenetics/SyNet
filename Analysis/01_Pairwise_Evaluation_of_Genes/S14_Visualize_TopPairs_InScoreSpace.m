clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/OScore/');
addpath('../_Utilities/');
ge_name = 'SyNet';

%% Load Top pairs
net_name = ['./Network_Files/' 'DSN_' ge_name '.mat'];
dsn_info = load(net_name, 'Pair_AUC', 'Gene_Name');
% dsn_info.Gene_Name = dsn_info.Gene_Name(1:500); dsn_info.Pair_AUC = dsn_info.Pair_AUC(1:500,1:500); %###
n_gene = numel(dsn_info.Gene_Name);
n_total = n_gene*(n_gene-1)/2;
Pair_Info = zeros(n_total, 14, 'single');
[Pair_Info(:,1), Pair_Info(:,2)] = find(triu(ones(n_gene), 1));
fprintf('In total [%d] gene pairs exist.\n', n_total);

%% Compute initial axes
Ind_AUC = dsn_info.Pair_AUC(1:n_gene+1:end)';
Pair_Info(:,3) = Ind_AUC(Pair_Info(:,1));
Pair_Info(:,4) = Ind_AUC(Pair_Info(:,2));
Pair_Info(:,5) = max(Pair_Info(:,3:4), [], 2); % Max AUC
pair_ind = sub2ind([n_gene n_gene], Pair_Info(:,1), Pair_Info(:,2));
Pair_Info(:,6) = dsn_info.Pair_AUC(pair_ind); % Combined AUC
Pair_Info(:,7) = Pair_Info(:,6)./Pair_Info(:,5); % Synergy
Pair_Info(:,8) = mean(Pair_Info(:,3:4), 2); % Mean AUC
dsn_info.Pair_AUC = [];
clear Ind_AUC

%% Add absolute correlation
GE_Path = getPath(ge_name);
fprintf('Loading [%s] expression data\n', GE_Path);
ge_info = load(GE_Path, 'Gene_Expression', 'Gene_Name');
% ge_info.Gene_Name = ge_info.Gene_Name(1:500); ge_info.Gene_Expression = ge_info.Gene_Expression(:,1:500); %###
if ~isequal(ge_info.Gene_Name, dsn_info.Gene_Name), error(); end
Abs_Crr = abs(corr(zscore(ge_info.Gene_Expression), 'Type', 'Spearman'));
Pair_Info(:,9) = Abs_Crr(pair_ind); % Absolute spearman correlation
clear pair_ind Abs_Crr
ge_info.Gene_Expression = [];

%% Normalizing scores
Axes_Name(7:9) = {'Synergy' 'Mean AUC' 'AbsCorr'};
Pair_Info(:,10) = oscore(Pair_Info(:,7)); % Synergy
Pair_Info(:,11) = oscore(Pair_Info(:,8)); % Mean AUC
Pair_Info(:,12) = oscore(Pair_Info(:,9)); % Absolute spearman correlation

%% Compute final fitness and sorting
for ai=10:12
    Pair_Info(:,13) = Pair_Info(:,13) + (1-Pair_Info(:,ai)).^2;
end
Pair_Info(:,13) = -sqrt(Pair_Info(:,13));
[~, sid] = sort(Pair_Info(:,13), 'Descend');
Pair_Info = Pair_Info(sid, :);
Pair_Info(:,14) = oscore(Pair_Info(:,13));

%% Plotting
close all
figure('Visible', 'off');
ind = [8 7];
n_top = 10000;
is_top = false(n_total,1); is_top(1:n_top) = 1;
plot(Pair_Info(~is_top, ind(1)), Pair_Info(~is_top, ind(2)), '.', 'Color', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerSize', 1);
hold on
plot(Pair_Info( is_top, ind(1)), Pair_Info( is_top, ind(2)), '.', 'Color', [0.3 0.8 0.2], 'MarkerEdgeColor', [0.3 0.8 0.2], 'MarkerSize', 4);
set(gca,'TickDir','out');

[hit_map, center] = hist3([Pair_Info(:,ind(1)) Pair_Info(:,ind(2))], 'Nbins', [50 50]);
cnt_h = contour(center{1}, center{2}, hit_map', 10, 'ShowText', 'on');
% cnt_h = contour(center{1}, center{2}, hit_map', 'LineStyle', '--', 'ShowText', 'on', 'LabelSpacing', inf, 'LevelList', logspace(1,7,7));
clr_map = gray(25); clr_map = clr_map(17:end,:);
colormap(clr_map);

set(gca, 'FontWeight', 'Bold');
xlabel(Axes_Name{ind(1)});
ylabel(Axes_Name{ind(2)});
title(['Score Space (' ge_name  ')']);
x_lim = [min(Pair_Info(:,ind(1)))*0.98 max(Pair_Info(:,ind(1)))*1.02];
y_lim = [min(Pair_Info(:,ind(2)))*0.98 max(Pair_Info(:,ind(2)))*1.02];
xlim(x_lim);
ylim(y_lim);

%% Saving the plot
sav_name = sprintf('./Plots/ScoreSpace_%s_TripleScores_%s-%s.png', ge_name, Axes_Name{ind(1)}, Axes_Name{ind(2)});
print(gcf, '-dpng', '-r300', sav_name);

%% Save the network
tsv_name = sprintf([net_name(1:end-4) '_TopPairs_N%d.tsv'], n_top);
fid = fopen(tsv_name, 'w');
fprintf(fid, 'Source\tTarget\tSynergy\tMean_AUC\tAbsCorr\tFitness\tWeight\n');
for pi=1:n_top
    fprintf(fid, '%s\t%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.5f\n', ge_info.Gene_Name{Pair_Info(pi,1)}, ge_info.Gene_Name{Pair_Info(pi,2)}, Pair_Info(pi,[7 8 9 13 14]));
end
fclose(fid);

