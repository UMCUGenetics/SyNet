function Pair_Info = getTopPairs(Target_Study, n_Pairs, Study_Index)


%% Initialization
dsn_path = '../01_Pairwise_Evaluation_of_Genes/Network_Files/';

%% Loading DS network
dsn_name = sprintf([dsn_path 'DSN_SyNetS%02d.mat'], Target_Study);
load(dsn_name, 'Net_Adj', 'anc_data');
if ~issymmetric(Net_Adj), error(); end
n_gene = size(Net_Adj, 1);
% [n_Study, n_Rep] = size(anc_data.cv_obj);

%% Sanity check
Used_Indice = unique(Study_Index(anc_data.cv_obj(Target_Study,1).iTr | anc_data.cv_obj(Target_Study,1).iTe));
if any(ismember(Used_Indice, Target_Study))
	error();
end

%% Sorting
[Pair_val, Pair_sid] = sort(Net_Adj(:), 'Descend');
Pair_sid = Pair_sid(1:n_Pairs*3);
Pair_val = Pair_val(1:n_Pairs*3);
Pair_Info = zeros(n_Pairs*3, 3);
[Pair_Info(:,1), Pair_Info(:,2)] = ind2sub([n_gene, n_gene], Pair_sid);
Pair_Info(:,3) = Net_Adj(Pair_sid);
if ~isequal(Pair_val, Pair_Info(:,3)), error(); end

%% Filtering
Pair_Info(Pair_Info(:,1)>Pair_Info(:,2), :) = [];
Pair_Info = Pair_Info(1:n_Pairs, :);
end