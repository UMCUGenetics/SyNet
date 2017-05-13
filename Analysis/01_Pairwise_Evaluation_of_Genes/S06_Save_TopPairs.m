clc;
clear;

%% Initialization
ge_name = 'SyNet';
net_path = './Network_Files/';

%% Load gene name
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
clear pair_ind Pair_AUC

%% Normalizing scores
ox = Pair_Info(:,6);
ox = ox-min(ox);
ox = ox/max(ox);

oy = Pair_Info(:,7);
oy = oy-min(oy);
oy = oy/max(oy);

Pair_Info(:,8) = -sqrt((1-ox).^2 + (1-oy).^2);

%% Sorting
[~, sid] = sort(Pair_Info(:,8), 'Descend');
Pair_Info = Pair_Info(sid, :);

%% Saving top pairs
n_study = 14;
PP_Info = Pair_Info(1:100000, :);
n_PP = size(PP_Info,1);
PP_PerStudy = zeros(n_PP, n_study+2, 'single');

% rid = randperm(n_total-n_PP, n_PP*5) + n_PP;
rid = randperm(n_total, n_PP*3);
NP_Info = Pair_Info(rid, :);
n_NP = size(NP_Info,1);
NP_PerStudy = zeros(n_NP, n_study+2, 'single');

%% Get per study score
PP_ID = arrayfun(@(i) sprintf('%d-%d', PP_Info(i,1:2)), 1:n_PP, 'UniformOutput', 0)';
PP_Map = containers.Map(PP_ID, 1:n_PP);
if numel(PP_ID) ~= PP_Map.Count, error(); end

NP_ID = arrayfun(@(i) sprintf('%d-%d', NP_Info(i,1:2)), 1:n_NP, 'UniformOutput', 0)';
NP_Map = containers.Map(NP_ID, 1:n_NP);
if numel(NP_ID) ~= NP_Map.Count, error(); end

pwr_path = ['./PWR_Files/' ge_name '/'];
pwr_lst = dir([pwr_path 'PWR_' ge_name '_*.mat']);
n_pwr = numel(pwr_lst);
fprintf('Found [%d] results.\n', n_pwr);
for ri=1:n_pwr
	res_name = [pwr_path pwr_lst(ri).name];
	fprintf('[%d/%d] Loading [%s]. ', ri, n_pwr, res_name);
	res_data = load(res_name, 'auc_pair');
	is_in = ismember(res_data.auc_pair(:,1:2), PP_Info(:,1:2), 'rows') | ismember(res_data.auc_pair(:,1:2), NP_Info(:,1:2), 'rows');
	if any(is_in)
		in_lst = find(is_in);
		fprintf('Found [%d] pairs.\n', numel(in_lst));
		for i=1:numel(in_lst)
			Item_ID = sprintf('%d-%d', res_data.auc_pair(in_lst(i),1:2));
			if PP_Map.isKey(Item_ID)
				pr_ind = PP_Map(Item_ID);
				PP_PerStudy(pr_ind,:) = res_data.auc_pair(in_lst(i),:);
			end
			if NP_Map.isKey(Item_ID)
				pr_ind = NP_Map(Item_ID);
				NP_PerStudy(pr_ind,:) = res_data.auc_pair(in_lst(i),:);
			end
		end
	else
		fprintf('\n');
	end
end
save(['./Top_Pairs/Top_' ge_name '.mat'], 'Ind_AUC', 'NP_Info', 'PP_Info', 'Gene_Name', 'PP_PerStudy', 'NP_PerStudy');
