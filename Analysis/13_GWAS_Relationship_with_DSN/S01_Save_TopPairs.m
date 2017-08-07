clc;
clear;

%% Initialization
addpath('../_Utilities/');
ge_name = 'SyNet';
dsn_path = '../01_Pairwise_Evaluation_of_Genes/Network_Files/';

%% Load gene name
net_name = [dsn_path 'DSN_' ge_name '.mat'];
fprintf('Loading [%s]\n', net_name);
load(net_name, 'Pair_AUC', 'Gene_Name');
n_gene = size(Pair_AUC, 1);
n_total = n_gene*(n_gene-1)/2;
fprintf('In total [%d] gene pairs exist.\n', n_total);

%% Identify pairs
fprintf('Loading pair list.\n');
Pair_Info = zeros(n_total, 10, 'single');
[Pair_Info(:,1), Pair_Info(:,2)] = find(triu(ones(n_gene), 1));

%% Generate Pair Matrix
Ind_AUC = Pair_AUC(1:n_gene+1:end)';
Pair_Info(:,3) = Ind_AUC(Pair_Info(:,1));
Pair_Info(:,4) = Ind_AUC(Pair_Info(:,2));
Pair_Info(:,5) = max(Pair_Info(:,3:4), [], 2);
pair_ind = sub2ind([n_gene n_gene], Pair_Info(:,1), Pair_Info(:,2));
Pair_Info(:,6) = Pair_AUC(pair_ind);
Pair_Info(:,7) = mean(Pair_Info(:,3:4), 2);
Pair_Info(:,8) = Pair_Info(:,6)./Pair_Info(:,5);

ge_path = getPath(ge_name);
ge_data = load(getPath(ge_name), 'Gene_Expression');
ax_crr = abs(corr(ge_data.Gene_Expression, 'Type', 'Spearman'));
ax_crr(1:size(ax_crr,1)+1:end) = 0;
Pair_Info(:,9) = ax_crr(pair_ind);
clear ge_data ax_crr pair_ind Pair_AUC

Pair_Info(:,10) = -sqrt((1-oscore(Pair_Info(:,7))).^2 + (1-oscore(Pair_Info(:,8))).^2 + (1-oscore(Pair_Info(:,9))).^2);

%% Sorting
[~, sid] = sort(Pair_Info(:,10), 'Descend');
Pair_Info = Pair_Info(sid, :);

%% Saving top pairs
PP_Info = Pair_Info(1:100000, :);
n_PP = size(PP_Info,1);

% rid = randperm(n_total-n_PP, n_PP*5) + n_PP;
rid = randperm(n_total, n_PP*3);
NP_Info = Pair_Info(rid, :);
n_NP = size(NP_Info,1);

save(['./Top_Pairs/Top_' ge_name '.mat'], 'Ind_AUC', 'NP_Info', 'PP_Info', 'Gene_Name');

%% //////////////// Functions
function lst = oscore(lst)
lst = lst - min(lst(:));
lst = lst ./ max(lst(:));
end