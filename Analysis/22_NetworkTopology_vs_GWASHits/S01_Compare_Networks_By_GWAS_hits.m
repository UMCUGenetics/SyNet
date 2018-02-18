clc;
clear
close all

%% Initialization
addpath('../_Utilities/');
% MAX_SyNet_Pairs = 03544;
MAX_SyNet_Pairs = 10000;
Ref_Name = 'AvgSynACr';
MAX_DISTANCE = 10000;
MAX_N_SNP = 53837; % pval = 0.00499986
% MAX_N_SNP = 10000; % pval = 0.000999939
% MAX_N_SNP = 5000;
% MAX_N_SNP = 3000;
% MAX_N_SNP = 1000;

%% Get reference gene set
net_opt.GE_Path = getPath('SyNet');
Ref_Data = load(net_opt.GE_Path, 'Gene_Name');
Reference_Gene_List = Ref_Data.Gene_Name;
n_ref_gene = numel(Reference_Gene_List);
Ref_GMap = containers.Map(Reference_Gene_List, 1:n_ref_gene);
clear Ref_Data

%% Load SyNet
dsn_name = ['../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet_' Ref_Name '.mat'];
fprintf('Loading [%s] network.\n', dsn_name);
DSN_info = load(dsn_name);
if ~isequal(DSN_info.Gene_Name, Reference_Gene_List), error(); end
SyNet_PairIndex = DSN_info.PP_Info(1:MAX_SyNet_Pairs,1:2);
SyNet_Gene_Index = unique(SyNet_PairIndex', 'Stable');
SyNet_Gene_Name = DSN_info.Gene_Name(SyNet_Gene_Index);
n_SyNet_Gene = numel(SyNet_Gene_Name);
ind = sub2ind([n_ref_gene, n_ref_gene], SyNet_PairIndex(:,1), SyNet_PairIndex(:,2));
SyNet_Adj = false(n_ref_gene);
SyNet_Adj(ind) = 1;
SyNet_Adj = max(SyNet_Adj, SyNet_Adj');
SyNet_Graph = graph(SyNet_Adj, DSN_info.Gene_Name);
clear SyNet_Adj DSN_info

%% Load GWAS
gwas_name = sprintf('../13_GWAS_Relationship_with_DSN/iCOGS_Hits/iCOGS_Hits_Genes_NTS%0.1fk_MD%0.1fk.tsv', MAX_N_SNP/1e3, MAX_DISTANCE/1e3);
fprintf('Loading GWAS data from [%s]\n', gwas_name);
fid = fopen(gwas_name, 'r');
% Id	-Log10(pval)	#Hit	#Hit/Size
GWAS_Info = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
fclose(fid);
if ~isequal(GWAS_Info{1}, Reference_Gene_List), error(); end
GWAS_Info = [GWAS_Info{2:end}];
Hot_Gene_List = find(GWAS_Info(:,1)>0);
n_hot = numel(Hot_Gene_List);
Cld_Gene_List = setdiff(1:n_ref_gene, Hot_Gene_List)';
Cld_Gene_List = Cld_Gene_List(randperm(numel(Cld_Gene_List), n_hot));
if numel(intersect(Hot_Gene_List, Cld_Gene_List)), error(); end

%% Degree
Index = 3;
switch Index
    case 1
        Hot_dist = SyNet_Graph.degree(Reference_Gene_List(Hot_Gene_List));
        Cld_dist = SyNet_Graph.degree(Reference_Gene_List(Cld_Gene_List));
    case 1.1
        Hot_dist = GWAS_Info(SyNet_Gene_Index, 1);
        Cld_dist = GWAS_Info(randperm(n_ref_gene, n_SyNet_Gene), 1);
    case 2
        Hot_dist = SyNet_Graph.distances(Reference_Gene_List(Hot_Gene_List), Reference_Gene_List(Hot_Gene_List));
        Cld_dist = SyNet_Graph.distances(Reference_Gene_List(Cld_Gene_List), Reference_Gene_List(Cld_Gene_List));
    case 3
        Hot_dist = zeros(n_SyNet_Gene, 1);
        Cld_dist = zeros(n_SyNet_Gene, 1);
        for gi=1:n_SyNet_Gene
            Nei_lst = SyNet_Graph.neighbors(SyNet_Gene_Name(gi));
            n_nei = numel(Nei_lst);
            Nei_Index = cellfun(@(item) Ref_GMap(item), Nei_lst);
            Hot_dist(gi) = sum(GWAS_Info(Nei_Index, 1));
            Cld_dist(gi) = sum(GWAS_Info(randperm(n_ref_gene, n_nei), 1));
        end
end

%% Filtering
Hot_dist = nonzeros(Hot_dist);
Cld_dist = nonzeros(Cld_dist);
Hot_dist(isinf(Hot_dist)) = [];
Cld_dist(isinf(Cld_dist)) = [];

%% Plotting
figure();
hold on
BoxPlotEx(Hot_dist(:), 'Positions', 1);
BoxPlotEx(Cld_dist(:), 'Positions', 2);
xlim([0 3]);


