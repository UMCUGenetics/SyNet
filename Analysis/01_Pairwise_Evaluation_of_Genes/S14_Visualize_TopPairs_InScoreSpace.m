clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/OScore/');
addpath('../_Utilities/');
ge_name = 'SyNet';

%% Load Top pairs
net_name = ['./Network_Files/' 'DSN_' ge_name '.mat'];
dsn_info = load(net_name, 'Pair_AUC', 'Gene_Name');
dsn_info.Gene_Name = dsn_info.Gene_Name(1:100); dsn_info.Pair_AUC = dsn_info.Pair_AUC(1:100,1:100); %###
n_gene = numel(dsn_info.Gene_Name);
n_total = n_gene*(n_gene-1)/2;
Pair_Info = zeros(n_total, 10, 'single');
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
Pair_Info(:,8) = mean(Pair_Info(:,3:4), 2); % Average AUC
dsn_info.Pair_AUC = [];

%% Add absolute correlation
GE_Path = getPath(ge_name);
fprintf('Loading [%s] expression data\n', GE_Path);
ge_info = load(GE_Path, 'Gene_Expression', 'Gene_Name');
ge_info.Gene_Name = ge_info.Gene_Name(1:100); ge_info.Gene_Expression = ge_info.Gene_Expression(:,1:100); %###
if ~isequal(ge_info.Gene_Name, dsn_info.Gene_Name), error(); end
Abs_Crr = corr(zscore(ge_info.Gene_Expression), 'Type', 'Spearman');
Pair_Info(:,9) = Abs_Crr(pair_ind); % Absolute spearman correlation
clear pair_ind Abs_Crr
ge_info.Gene_Expression = [];

%% Normalizing scores
Pair_Info(:,7) = oscore(Pair_Info(:,7));
Pair_Info(:,8) = oscore(Pair_Info(:,8));
Pair_Info(:,9) = oscore(Pair_Info(:,9));

%% Compute final fitness and sorting
for ai=[7 8 9]
    Pair_Info(:,10) = (1-Pair_Info(:,ai)).^2;
end
Pair_Info(:,10) = -sqrt(Pair_Info(:,10));
[~, sid] = sort(Pair_Info(:,10), 'Descend');
Pair_Info = Pair_Info(sid, :);

%% Plotting
close all
figure();
plot(NP_Info(:,ind(1)), NP_Info(:,ind(2)), '.', 'Color', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerSize', 1);
hold on
plot(PP_Info(:,ind(1)), PP_Info(:,ind(2)), '.', 'Color', [0.3 0.8 0.2], 'MarkerEdgeColor', [0.3 0.8 0.2], 'MarkerSize', 4);

[hit_map, center] = hist3([NP_Info(:,ind(1)), NP_Info(:,ind(2))], 'Nbins', [50 50]);
cnt_h = contour(center{1}, center{2}, hit_map', 10, 'ShowText', 'off');
clr_map = gray(25); clr_map = clr_map(17:end,:);
colormap(clr_map);

set(gca, 'FontWeight', 'Bold');
xlabel('Max AUC');
ylabel('Synergy');
title(['Score Space (' ge_name  ')']);
x_lim = [min(NP_Info(:,ind(1)))*0.98 max(NP_Info(:,ind(1)))*1.02];
y_lim = [min(NP_Info(:,ind(2)))*0.98 max(NP_Info(:,ind(2)))*1.02];
xlim(x_lim);
ylim(y_lim);

%% Saving the plot
sav_name = ['./Plots/ScoreSpace_' ge_name '_TripleScores.png'];
print(gcf, '-dpng', '-r300', sav_name);

