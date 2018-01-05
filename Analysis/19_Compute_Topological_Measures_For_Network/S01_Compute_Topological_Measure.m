function S01_Compute_Topological_Measure(net_name, TM_Name, RemoveGSet)
% clc;
%{
for ti in DirectConnection ShortestPath Degree PageRank-FB0.65 PageRank-FB0.75 PageRank-FB0.85 PageRank-FB0.95 Closeness Betweenness Eigenvector; do
for ni in HumanInt HBBrain HBKidney HBOvary HBUterus AbsCorr STRING IntAct BioPlex BioGRID HBLympNode HBEpith HBGland HBOvary; do
PARAM=\'$ni\',\'$ti\';
sbatch --job-name=CTM-$PARAM --output=Logs/CTM-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=10GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S01_Compute_Topological_Measure "$PARAM"; 
done; 
read -p "`date`: $PARAM. Press a key" -t 1800
done
%}

%% Inialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
if ~exist('RemoveGSet', 'var'), RemoveGSet={}; end
if ismac
    net_name = 'STRING'; % 'I2D' 'STRING' 'HPRD' 'HBEpith','HBGland'
    TM_Name = 'DirectConnection';
end
Ref_Name = 'AvgSynACr';
MAX_SyNet_Pair = 20000;
Shuff_Method = 'LnkShuff';

%% Load Pair index
PI_Name = sprintf('./Topological_Data/PairInfo-%s_%s_MP%06d.mat', Shuff_Method, Ref_Name, MAX_SyNet_Pair);
fprintf('Loading pair index info [%s]\n', PI_Name);
Pair_Data = load(PI_Name);
Ref_GeneName = Pair_Data.SyNet_info.Gene_Name;
Pair_Info = [
    Pair_Data.PP_Info(:, 1:2)
    Pair_Data.NP_Info
    ];
n_RefGeneName = numel(Ref_GeneName);
clear Pair_Data

%% Load network
fprintf('Loading [%s] as existing network.\n', net_name);
net_opt.GE_Path = getPath('SyNet');
net_opt.PreferredGenes = Ref_GeneName;
net_opt.MAX_N_PAIR = 50000;
net_info = LoadNetworkAdj(net_name, net_opt);
NET_N_PAIR = net_info.N_PAIR;
if min(net_info.Net_Adj(:))<0, error('Not implemented for negative links'); end

%% Filter network
if ~isempty(RemoveGSet)
    is_in = ismember(net_info.Gene_Name, RemoveGSet);
    fprintf('[%d] genes are given: [%d] genes from total of [%d] genes are filtered.\n', nueml(RemoveGSet), sum(is_in), numel(is_in));
    net_info.Gene_Name(is_in) = [];
    net_info.Net_Adj(is_in, :) = [];
    net_info.Net_Adj(:, is_in) = [];
end

%% Reorder Adjacency matrix
fprintf('Reordering adjacency matrix [%s] according to [%s] genes.\n', net_name, Ref_Name);
GMap = containers.Map(net_info.Gene_Name, 1:numel(net_info.Gene_Name));
Ref2Net = zeros(n_RefGeneName, 2);
for gi=1:n_RefGeneName
    if GMap.isKey(Ref_GeneName{gi})
        Ref2Net(gi,:) = [gi GMap(Ref_GeneName{gi})];
    end
end
Ref2Net(Ref2Net(:,1)==0, :) = [];
if ~isequal(Ref_GeneName(Ref2Net(:,1)), net_info.Gene_Name(Ref2Net(:,2))), error(); end
Net_Adj = zeros(n_RefGeneName, 'single');
Net_Adj(Ref2Net(:,1), Ref2Net(:,1)) = net_info.Net_Adj(Ref2Net(:,2), Ref2Net(:,2));
clear net_info
fprintf('[%d] genes and [%d] links are left in the network.\n', size(Net_Adj,1), numel(nonzeros(triu(Net_Adj))));

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
Output_results = struct();
Output_results.Ref_GeneName = Ref_GeneName;
Output_results.Pair_Info = Pair_Info;
Output_results.NET_N_PAIR = NET_N_PAIR;
Output_results.net_opt = net_opt;
fprintf('Computing [%s] from network [%s] ...\n', TM_Name, net_name);
Net_Graph = graph(Net_Adj~=0, 'OmitSelfLoops');
Pair_Index = sub2ind([n_RefGeneName n_RefGeneName], Pair_Info(:,1), Pair_Info(:,2));
switch TM_Name
    case 'DirectConnection'
        Output_results.Pair_AvgScore = double(Net_Adj(Pair_Index)>0);
        Output_results.Pair_DifScore = Output_results.Pair_AvgScore;
    case 'ShortestPath'
        Dist_Mat = Net_Graph.distances('Method', 'unweighted');
        Output_results.Pair_AvgScore = Dist_Mat(Pair_Index);
        Output_results.Pair_DifScore = Output_results.Pair_AvgScore;
    case 'Degree'
        Output_results.Gene_Degree = centrality(Net_Graph, 'degree');
        [Output_results.Pair_AvgScore, Output_results.Pair_DifScore] = ConvertGeneToPairs(Output_results.Gene_Degree, Pair_Info(:,1:2));
    case {'PageRank-FB0.65' 'PageRank-FB0.75' 'PageRank-FB0.85' 'PageRank-FB0.95'}
        Output_results.FollowProbability = str2double(TM_Name(end-3:end));
        Output_results.Gene_PageRank = centrality(Net_Graph, 'pagerank', 'FollowProbability', Output_results.FollowProbability);
        [Output_results.Pair_AvgScore, Output_results.Pair_DifScore] = ConvertGeneToPairs(Output_results.Gene_PageRank, Pair_Info(:,1:2));
    case 'Closeness'
        Output_results.Gene_Closeness = centrality(Net_Graph, 'closeness');
        [Output_results.Pair_AvgScore, Output_results.Pair_DifScore] = ConvertGeneToPairs(Output_results.Gene_Closeness, Pair_Info(:,1:2));
    case 'Betweenness'
        Output_results.Gene_Betweenness = centrality(Net_Graph, 'betweenness');
        [Output_results.Pair_AvgScore, Output_results.Pair_DifScore] = ConvertGeneToPairs(Output_results.Gene_Betweenness, Pair_Info(:,1:2));
    case 'Eigenvector'
        Output_results.Gene_Eigenvector = centrality(Net_Graph, 'eigenvector');
        [Output_results.Pair_AvgScore, Output_results.Pair_DifScore] = ConvertGeneToPairs(Output_results.Gene_Eigenvector, Pair_Info(:,1:2));
end

%% Save the results
sav_name = sprintf('./Topological_Results/TM-%s_%s_%s_NP%06d_MSP%06d_%s.mat', Shuff_Method, Ref_Name, net_name, net_opt.MAX_N_PAIR, MAX_SyNet_Pair, TM_Name);
fprintf('Saving the results in [%s]\n', sav_name);
save(sav_name, '-struct', 'Output_results');
end

%% Functions %%%%%%%%%%%%%%%%%%%%%%
function [Avg_Score, Diff_Score] = ConvertGeneToPairs(Gene_Score, Pair_Indices)
Pair_Score = sort([Gene_Score(Pair_Indices(:,1)) Gene_Score(Pair_Indices(:,2))], 2, 'Descend');
Avg_Score = mean(Pair_Score, 2);
Diff_Score = Pair_Score(:,1) - Pair_Score(:,2);
end

