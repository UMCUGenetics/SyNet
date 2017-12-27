clc;
clear;

%% Initialization
outout_path = './Cyto_Input/';
sn_name = 'SyNet';
MAX_TOPPAIR = 3544;

%% Load SyNet pairs
dsn_name = '../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet_AvgSynACr.mat';
fprintf('Loading [%s] network.\n', dsn_name);
DSN_info = load(dsn_name);
SyNet_PairIndex = DSN_info.PP_Info(1:MAX_TOPPAIR,1:2);
SyNet_PairWeight = DSN_info.PP_Info(1:MAX_TOPPAIR,15);
All_GeneName = DSN_info.Gene_Name;

%% Output nodes
Top_NodeName = All_GeneName(unique(SyNet_PairIndex', 'stable'));
n_node = numel(Top_NodeName);
node_fname = sprintf([outout_path sn_name '_Top%05d_Node.tsv'], MAX_TOPPAIR);
fprintf('Writing [%d] node output to [%s].\n', n_node, node_fname);
fid = fopen(node_fname, 'w');
fprintf(fid, 'Id\tLabel\tWeight\n');
for ni=1:n_node
    fprintf(fid, '%s\t%s\t%d\n', Top_NodeName{ni}, Top_NodeName{ni}, n_node-ni+1);
end
fclose(fid);

%% Output edges
n_edge = size(SyNet_PairIndex, 1);
edge_fname = sprintf([outout_path sn_name '_Top%05d_Edge.tsv'], MAX_TOPPAIR);
fprintf('Writing [%d] edge output to [%s].\n', n_edge, edge_fname);
fid = fopen(edge_fname, 'w');
fprintf(fid, 'Source\tTarget\tType\tWeight\n');
for ei=1:n_edge
    fprintf(fid, '%s\t%s\tUndirected\t%f\n', All_GeneName{SyNet_PairIndex(ei,1)}, All_GeneName{SyNet_PairIndex(ei,2)}, SyNet_PairWeight(ei));
end
fclose(fid);
