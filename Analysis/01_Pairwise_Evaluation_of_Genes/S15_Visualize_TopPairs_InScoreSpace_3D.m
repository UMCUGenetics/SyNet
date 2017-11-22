clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/OScore/');
addpath('../_Utilities/');
ge_name = 'SyNet';
n_top = 10000;

%% Load Top pairs
net_name = ['./Network_Files/' 'DSN_' ge_name '.mat'];
dsn_info = load(net_name, 'Pair_AUC', 'Gene_Name');
% dsn_info.Gene_Name = dsn_info.Gene_Name(1:500); dsn_info.Pair_AUC = dsn_info.Pair_AUC(1:500,1:500); %###
n_gene = numel(dsn_info.Gene_Name);
n_total = n_gene*(n_gene-1)/2;
Pair_Info = zeros(n_total, 15, 'single');
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
% clear Ind_AUC

%% Add absolute correlation
GE_Path = getPath(ge_name);
fprintf('Loading [%s] expression data\n', GE_Path);
ge_info = load(GE_Path, 'Gene_Expression', 'Gene_Name');
% ge_info.Gene_Name = ge_info.Gene_Name(1:500); ge_info.Gene_Expression = ge_info.Gene_Expression(:,1:500); %###
if ~isequal(ge_info.Gene_Name, dsn_info.Gene_Name), error(); end
Corr_mat = corr(zscore(ge_info.Gene_Expression), 'Type', 'Spearman');
Pair_Info(:,9) = Corr_mat(pair_ind);
Pair_Info(:,10) = abs(Pair_Info(:,9)); % Absolute spearman correlation
clear pair_ind Corr_mat
ge_info.Gene_Expression = [];

%% Normalizing scores
Axes_Name([7 8 10]) = {'Synergy' 'Mean AUC' 'AbsCorr'};
Pair_Info(:,11) = oscore(Pair_Info(:,7)); % Synergy
Pair_Info(:,12) = oscore(Pair_Info(:,8)); % Mean AUC
Pair_Info(:,13) = oscore(Pair_Info(:,10)); % Absolute spearman correlation

%% Compute final fitness and sorting
for ai=11:13
    Pair_Info(:,14) = Pair_Info(:,14) + (1-Pair_Info(:,ai)).^2;
end
Pair_Info(:,14) = -sqrt(Pair_Info(:,14));
[~, sid] = sort(Pair_Info(:,14), 'Descend');
Pair_Info = Pair_Info(sid, :);
Pair_Info(:,15) = oscore(Pair_Info(:,14));

%% Plotting
close all
figure('Visible', 'off');
ind = [8 7 10];
is_top = false(n_total,1); is_top(1:n_top) = 1;
% rest_h = scatter3(Pair_Info(~is_top, ind(1)), Pair_Info(~is_top, ind(2)), Pair_Info(~is_top, ind(3)), 0.5, [0.8 0.8 0.8]);
% rest_h.MarkerEdgeAlpha = 0.2;
plot3(Pair_Info(~is_top, ind(1)), Pair_Info(~is_top, ind(2)), Pair_Info(~is_top, ind(3)), '.', 'Color', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerSize', 1);
% set(gca, 'YGrid', 'off', 'XGrid', 'off', 'ZGrid', 'off');
hold on
plot3(Pair_Info( is_top, ind(1)), Pair_Info( is_top, ind(2)), Pair_Info( is_top, ind(3)), '.', 'Color', [0.3 0.8 0.2], 'MarkerEdgeColor', [0.3 0.8 0.2], 'MarkerSize', 3);
set(gca,'TickDir','out', 'FontWeight', 'Bold', 'YGrid', 'on', 'XGrid', 'on', 'ZGrid', 'on', 'GridAlpha', 0.1, ...
    'XTick', 0.52:0.02:0.64, 'YTick', 0.85:0.05:1.15, 'ZTick', 0.1:0.2:1.0);

xlabel(Axes_Name{ind(1)});
ylabel(Axes_Name{ind(2)});
zlabel(Axes_Name{ind(3)});
title(['Score Space (' ge_name  ')']);
x_lim = [0.50 0.66]; %[min(Pair_Info(:,ind(1)))*0.98 max(Pair_Info(:,ind(1)))*1.02];
y_lim = [0.83 1.16]; %[min(Pair_Info(:,ind(2)))*0.98 max(Pair_Info(:,ind(2)))*1.02];
z_lim = [0 0.95]; %[min(Pair_Info(:,ind(3)))*0.98 max(Pair_Info(:,ind(3)))*1.02];
xlim(x_lim);
ylim(y_lim);
zlim(z_lim);
view(45,30);

%% Saving the plot
sav_name = sprintf('./Plots/ScoreSpace_%s_TripleScores_In3D_%s-%s-%s_3D.png', ge_name, Axes_Name{ind(1)}, Axes_Name{ind(2)}, Axes_Name{ind(3)});
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'portrait', 'PaperPositionMode','auto', 'PaperSize', [6 4], 'PaperPosition', [0 0 6 4]);
print(gcf, '-dpng', '-r300', sav_name);
