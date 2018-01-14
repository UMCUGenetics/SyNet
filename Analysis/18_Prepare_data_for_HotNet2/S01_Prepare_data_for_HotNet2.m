clc;
clear

%% Initialization
addpath('../17_Frequency_of_SyNetLinks_in_Biological_Netwrok/');
addpath('../_Utilities/');
MAX_DISTANCE = 10000;
MAX_SyNet_Pairs = 50000;
% MAX_N_SNP = 504752;
MAX_N_SNP = 53837; % pval = 0.00499986
% MAX_N_SNP = 10962; % pval = 0.000999939
% MAX_N_SNP = 10000; % pval = 0.000999939
% MAX_N_SNP = 1104; % pval = 0.0000999498
% MAX_N_SNP = 1000;

%% Load SyNet mat file
dsn_name = ['../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet_AvgSynACr.mat'];
fprintf('Loading [%s] network.\n', dsn_name);
DSN_info = load(dsn_name);
n_gene = numel(DSN_info.Gene_Name);

%% Load GWAS hits over Cohort
fid = fopen(sprintf('../13_GWAS_Relationship_with_DSN/iCOGS_Hits/iCOGS_Hits_Genes_NTS%0.1fk_MD%0.1fk.tsv', MAX_N_SNP/1e3, MAX_DISTANCE/1e3), 'r');
% Id	-Log10(pval)	#Hit	#Hit/Size
GWAS_Info = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
fclose(fid);
if ~isequal(GWAS_Info{1}, DSN_info.Gene_Name), error(); end

%% Output makeHeatFile.py pvalues --heat_file
GWAS_Info{2} = GWAS_Info{2} / max(GWAS_Info{2});
fid = fopen(sprintf('./HotNet2_Input_Files/iCOGS_NTS%0.1fk_MD%0.1fk_LogPval.tsv', MAX_N_SNP/1e3, MAX_DISTANCE/1e3), 'w');
for gi=1:n_gene
    fprintf(fid, '%s\t%0.5f\n', GWAS_Info{1}{gi}, GWAS_Info{2}(gi));
end
fclose(fid);

%% Output makeHeatFile.py # Hits notmalized by size --heat_file
GWAS_Info{4} = GWAS_Info{4} / max(GWAS_Info{4});
fid = fopen(sprintf('./HotNet2_Input_Files/iCOGS_NTS%0.1fk_MD%0.1fk_NHitSize.tsv', MAX_N_SNP/1e3, MAX_DISTANCE/1e3), 'w');
for gi=1:n_gene
    fprintf(fid, '%s\t%0.5f\n', GWAS_Info{1}{gi}, GWAS_Info{4}(gi));
end
fclose(fid);

%% Output SyNet
OutputHotNetwork('SyNet', DSN_info.PP_Info(1:MAX_SyNet_Pairs,1:2), DSN_info.Gene_Name);
% OutputHotNetwork('SyNet', DSN_info.PP_Info( 1:20000,1:2), DSN_info.Gene_Name);
% OutputHotNetwork('SyNet', DSN_info.PP_Info( 1:50000,1:2), DSN_info.Gene_Name);
% OutputHotNetwork('SyNet', DSN_info.PP_Info(1:100000,1:2), DSN_info.Gene_Name);

%% Output Shuffled SyNet
rind = randperm(n_gene);
OutputHotNetwork('SyNet-Shuff', DSN_info.PP_Info(1:MAX_SyNet_Pairs,1:2), DSN_info.Gene_Name(rind));
% OutputHotNetwork('SyNet-Shuff', DSN_info.PP_Info( 1:20000,1:2), DSN_info.Gene_Name(rind));
% OutputHotNetwork('SyNet-Shuff', DSN_info.PP_Info( 1:50000,1:2), DSN_info.Gene_Name(rind));
% OutputHotNetwork('SyNet-Shuff', DSN_info.PP_Info(1:100000,1:2), DSN_info.Gene_Name(rind));

%% Output STRING
% net_opt.GE_Path = getPath('SyNet');
% ge_data = load(net_opt.GE_Path, 'Gene_Name');
% net_opt.PreferredGenes = ge_data.Gene_Name;
% net_opt.MAX_N_PAIR = 150000;
% net_info = LoadNetworkPairs('STRING', net_opt);
% OutputHotNetwork('STRING', net_info.Pair_Index( 1:10000,1:2), net_info.Gene_Name);
% OutputHotNetwork('STRING', net_info.Pair_Index( 1:20000,1:2), net_info.Gene_Name);
% OutputHotNetwork('STRING', net_info.Pair_Index( 1:50000,1:2), net_info.Gene_Name);
% OutputHotNetwork('STRING', net_info.Pair_Index(1:100000,1:2), net_info.Gene_Name);

%% Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function OutputHotNetwork(Net_Name, Pair_List, Gene_Name)
n_pair = size(Pair_List, 1);
n_gene = numel(Gene_Name);

%% Output edge list: makeNetworkFiles.py --edgelist_file
fid = fopen(sprintf('./HotNet2_Input_Files/%s_NP%06d_edgelist.tsv', Net_Name, n_pair), 'w');
for pi=1:n_pair
    fprintf(fid, '%d\t%d\n', Pair_List(pi,1:2));
end
fclose(fid);

%% Output Index to Gene: makeNetworkFiles.py --gene_index_file
fid = fopen(sprintf('./HotNet2_Input_Files/%s_geneindex.tsv', Net_Name), 'w');
for gi=1:n_gene
    fprintf(fid, '%d\t%s\n', gi, Gene_Name{gi});
end
fclose(fid);
end

function net_info = LoadNetworkPairs(net_name, net_opt)
switch net_name
    case {'STRING','HPRD','I2D','HBEpith','HBGland'}
        net_info.net_path = getPath(net_name);
        fid = fopen(net_info.net_path, 'r');
        Header_lst = regexp(fgetl(fid), '\t', 'split');
        if numel(Header_lst)==2
            fprintf('No link weight exists. Selection of [%d] links from begining of file.\n', net_opt.MAX_N_PAIR*5);
            net_cell = textscan(fid,   '%s%s', net_opt.MAX_N_PAIR, 'Delimiter', '\t', 'ReturnOnError', 0);
        else
            fprintf('Link weights are found. Selecting of [%d] links from top weighted interactions.\n', net_opt.MAX_N_PAIR*5);
            net_cell = textscan(fid, '%s%s%d', net_opt.MAX_N_PAIR, 'Delimiter', '\t', 'ReturnOnError', 0);
            if ~issorted(net_cell{3}, 'Descend'), error('Given network file is not sorted.'); end
        end
        fclose(fid);
        is_val = ismember(net_cell{1}, net_opt.PreferredGenes) & ismember(net_cell{2}, net_opt.PreferredGenes);
        net_cell = {net_cell{1}(is_val) net_cell{2}(is_val)};
    otherwise
        error('Unknown method.');
end

%% Convering gene names fo indices
n_lnk = numel(net_cell{1});
n_gene = numel(net_opt.PreferredGenes);
GMap = containers.Map(net_opt.PreferredGenes, 1:n_gene);
net_info.Pair_Index = zeros(n_lnk, 2, 'uint32');
for li=1:n_lnk
    net_info.Pair_Index(li,:) = [GMap(net_cell{1}{li}) GMap(net_cell{2}{li})];
end
net_info.Gene_Name = net_opt.PreferredGenes;
end