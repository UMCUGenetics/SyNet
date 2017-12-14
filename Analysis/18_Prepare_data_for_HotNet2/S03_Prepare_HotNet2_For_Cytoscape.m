clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/JSONlab/jsonlab-1.5/');
viz_file = './HotNet2_results/synet_np010000-icogs_10k/viz-data.json';

%% Load SyNet pairs
dsn_name = '../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat';
fprintf('Loading [%s] network.\n', dsn_name);
DSN_info = load(dsn_name);
SyNet_PairIndex = DSN_info.PP_Info(1:10000,1:2);
SyNet_GeneName = DSN_info.Gene_Name;
n_gene = numel(SyNet_GeneName);
ind = sub2ind([n_gene, n_gene], SyNet_PairIndex(:,1), SyNet_PairIndex(:,2));
SyNet_Adj = false(n_gene);
SyNet_Adj(ind) = 1;
SyNet_Adj = max(SyNet_Adj, SyNet_Adj');
SyNet_Graph = graph(SyNet_Adj, SyNet_GeneName);
clear SyNet_Adj

%% Loading json
viz_data = loadjson(viz_file);
opt_delta_str = viz_data.params.auto_delta;
Opt_Subnet_List = viz_data.subnetworks.(sprintf('x0x30__0x2E_%s', opt_delta_str(3:end)));

%% Output subnetworks
for si=1:numel(Opt_Subnet_List)
    sn_name = sprintf('SN_%03d_%s', si, opt_delta_str);
    fprintf('Processing [%s]: \n', sn_name);
    [HotNet_Gene_lst, HotNet_Pair_lst] = Output_SubNet(sn_name, Opt_Subnet_List{si}, viz_data.gene2heat);
    HotNet_n_pair = size(HotNet_Pair_lst, 1);
    HotNet_n_gene = numel(HotNet_Gene_lst);
    HotNet_Adj = false(HotNet_n_gene);
    for pi=1:HotNet_n_pair
        src_ind = find(strcmp(HotNet_Pair_lst{pi,1}, HotNet_Gene_lst));
        tar_ind = find(strcmp(HotNet_Pair_lst{pi,2}, HotNet_Gene_lst));
        if isempty(src_ind) || isempty(tar_ind), error(); end
        if src_ind==tar_ind, error(); end
        HotNet_Adj(src_ind, tar_ind) = 1;
        HotNet_Adj(tar_ind, src_ind) = 1;
    end
    HotNet_Graph = graph(HotNet_Adj, HotNet_Gene_lst);
    
    Extra_lnk = {};
    for gi=1:HotNet_n_gene
        HotNet_DiGraph = HotNet_Graph.shortestpathtree(HotNet_Gene_lst{gi}, HotNet_Gene_lst);
        if isempty(HotNet_DiGraph.Edges)
            fprintf('Gene [%s] is missing connection.\n', HotNet_Gene_lst{gi});
            HotNet_OtherGenes = setdiff(HotNet_Gene_lst, HotNet_Gene_lst{gi});
            [SyNet_DiGraph, nHops_To_SyNet] = SyNet_Graph.shortestpathtree(HotNet_Gene_lst{gi}, HotNet_OtherGenes);
            %del_node = find(centrality(SyNet_DiGraph, 'indegree') == 0 & centrality(SyNet_DiGraph, 'outdegree') == 0);
            %SyNet_DiGraph = SyNet_DiGraph.rmnode(del_node);
            MAX_N_HOPS = min(nHops_To_SyNet);
            for gj=1:numel(HotNet_OtherGenes)
                [hop_lst, n_hop] = SyNet_Graph.shortestpath(HotNet_Gene_lst{gi}, HotNet_OtherGenes{gj});
                if n_hop <= MAX_N_HOPS
                    for hi=1:n_hop-1
                        Extra_lnk = [Extra_lnk; sort(hop_lst(li:li+1))];
                    end
                end
            end
        end
    end
    Extra_lnk = unique(Extra_lnk, 'rows');
    
    edge_fname = ['./Cyto_Input/' sn_name '_Edge.tsv'];
    fid = fopen(edge_fname, 'a');
    for li=1:size(Extra_lnk,1)
        fprintf(fid, '%s\t%s\tUndirected\tDisconnected\n', Extra_lnk{li,1}, Extra_lnk{li,2});
    end
    fclose(fid);
end

function [Gene_lst, Pair_lst] = Output_SubNet(sn_name, sn_data, gene_heat)
outout_path = './Cyto_Input/';

%% Output nodes
node_fname = [outout_path sn_name '_Node.tsv'];
n_node = numel(sn_data.nodes);
Gene_lst = cell(n_node, 1);
fprintf('Writing [%d] node output to [%s].\n', n_node, node_fname);
fid = fopen(node_fname, 'w');
fprintf(fid, 'Id\tLabel\tGeneHeat\tHotNet2_HeatValue\n');
for ni=1:n_node
    Gene_lst{ni} = sn_data.nodes{ni}.name;
    gene_name = strrep(Gene_lst{ni}, '-', '_0x2D_');
    fprintf(fid, '%s\t%s\t%f\t%f\n', sn_data.nodes{ni}.name, sn_data.nodes{ni}.name, gene_heat.(gene_name), sn_data.nodes{ni}.value);
end
fclose(fid);

%% Output edges
edge_fname = [outout_path sn_name '_Edge.tsv'];
n_edge = numel(sn_data.edges);
Pair_lst = cell(n_edge, 2);
fprintf('Writing [%d] edge output to [%s].\n', n_edge, edge_fname);
fid = fopen(edge_fname, 'w');
fprintf(fid, 'Source\tTarget\tType\tCategories\n');
for ei=1:n_edge
    if numel(sn_data.edges{ei}.categories)~=1, error(); end
    fprintf(fid, '%s\t%s\tUndirected\t%s\n', sn_data.edges{ei}.source, sn_data.edges{ei}.target, sn_data.edges{ei}.categories{1});
    Pair_lst(ei,:) = {sn_data.edges{ei}.source sn_data.edges{ei}.target};
end
fclose(fid);
end