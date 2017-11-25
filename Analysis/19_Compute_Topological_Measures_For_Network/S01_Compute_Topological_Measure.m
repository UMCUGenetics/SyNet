clc;
clear;

%% Inialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
net_name = 'STRING';
n_pair = 2000;
tm_name = 'ShortestPath';

%% Load SyNet
SyNet_info = load('../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat', 'PP_Info', 'NP_Info', 'Gene_Name');
Pair_Info = [
    SyNet_info.PP_Info(1:n_pair,:)  ones(n_pair, 1);
    SyNet_info.NP_Info(1:n_pair,:) -ones(n_pair, 1);
    ];
Ref_GeneName = SyNet_info.Gene_Name;
n_RefGeneName = numel(Ref_GeneName);
clear SyNet_info

%% Load network
fprintf('Loading [%s] network.\n', net_name);
net_opt.PreferredGenes = Ref_GeneName;
net_opt.MAX_N_PAIR = n_pair;
net_info = LoadNetworkAdj(net_name, net_opt);
NET_N_PAIR = net_info.N_PAIR;
if min(net_info.Net_Adj(:))<0, error('Not implemented for negative links'); end

%% Reorder Adjacency matrix
GMap = containers.Map(net_info.Gene_Name, 1:numel(net_info.Gene_Name));
Ref2Net = zeros(n_RefGeneName, 2);
for gi=1:n_RefGeneName
    if GMap.isKey(Ref_GeneName{gi})
        Ref2Net(gi,:) = [gi GMap(Ref_GeneName{gi})];
    end
end
Ref2Net(Ref2Net(:,1)==0, :) = [];
Net_Adj = zeros(n_RefGeneName, 'single');
Net_Adj(Ref2Net(:,1), Ref2Net(:,1)) = net_info.Net_Adj(Ref2Net(:,2), Ref2Net(:,2));
clear net_info

%{
%%Test example
Net_Adj = zeros(6);
Ref_Adj = zeros(6);
for i=1:6
    for j=1:6
        Net_Adj(i,j) = i*10+j;
    end
end
i = [1 3 5];
r = [3 1 6];
Ref_Adj(r,r) = Net_Adj(i,i);
%}

%% Compute topological measiure
Output_results.Ref_GeneName = Ref_GeneName;
Output_results.Pair_Score = zeros(n_pair, 1);
fprintf('Computing [%s] from network [%s] ...\n', tm_name, net_name);
Net_Graph = graph(Net_Adj~=0);
Pair_Index = sub2ind([n_RefGeneName n_RefGeneName], Pair_Info(:,1), Pair_Info(:,2));
switch tm_name
    case 'Degree'
        Output_results.Gene_Degree = Net_Graph.degree;
        Pair_Degree = [Output_results.Gene_Degree(Pair_Info(:,1)) Output_results.Gene_Degree(Pair_Info(:,2))];
        Output_results.Pair_Score = mean(Pair_Degree, 2);
    case 'ShortestPath'
        Dist_Mat = Net_Graph.distances('Method', 'unweighted');
        Output_results.Pair_Score = Dist_Mat(Pair_Index);
end

%% Save the results
sav_name = sprintf('./Topological_Results/TM_%s_NP%06d_%s.mat', net_name, NET_N_PAIR, tm_name);
fprintf('Saving the results in [%s]\n', sav_name);
save(sav_name, '-struct', 'Output_results');



