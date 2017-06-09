clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
ge_name = 'SyNet';

%% Load Top pairs
ind = [6 7];
Pair_Info = collectTopPairs(ge_name, ind);
PP_Info = Pair_Info(1:10000,:);

% n_total = size(Pair_Info, 1);
% rid = randperm(n_total, n_PP*3);
NP_Info = Pair_Info;
clear Pair_Info

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
sav_name = ['./Plots/ScoreSpace_' ge_name '_CLR-CLS.png'];
print(gcf, '-dpng', '-r300', sav_name);

function Pair_Info = collectTopPairs(ge_name, ind)
%% Load gene name
net_path = './Network_Files/';
net_name = [net_path 'DSN_' ge_name '.mat'];
fprintf('Loading [%s]\n', net_name);
load(net_name, 'Pair_AUC', 'Gene_Name');
n_gene = size(Pair_AUC, 1);
n_total = n_gene*(n_gene-1)/2;
fprintf('In total [%d] gene pairs exist.\n', n_total);

%% Identify pairs
fprintf('Loading pair list.\n');
Pair_Info = zeros(n_total, 8, 'single');
[Pair_Info(:,1), Pair_Info(:,2)] = find(triu(ones(n_gene), 1));

%% Generate Pair Matrix
Ind_AUC = Pair_AUC(1:n_gene+1:end)';
Pair_Info(:,3) = Ind_AUC(Pair_Info(:,1));
Pair_Info(:,4) = Ind_AUC(Pair_Info(:,2));
Pair_Info(:,5) = max(Pair_Info(:,3:4), [], 2);
pair_ind = sub2ind([n_gene n_gene], Pair_Info(:,1), Pair_Info(:,2));
Pair_Info(:,6) = Pair_AUC(pair_ind);
Pair_Info(:,7) = Pair_Info(:,6)./Pair_Info(:,5);
Pair_Info(:,9) = mean(Pair_Info(:,3:4), 2);
clear pair_ind Pair_AUC

%% Normalizing scores
ox = Pair_Info(:,ind(1));
ox = ox-min(ox);
ox = ox/max(ox);

oy = Pair_Info(:,ind(2));
oy = oy-min(oy);
oy = oy/max(oy);

Pair_Info(:,8) = -sqrt((1-ox).^2 + (1-oy).^2);

%% Sorting
[~, sid] = sort(Pair_Info(:,8), 'Descend');
Pair_Info = Pair_Info(sid, :);
end
