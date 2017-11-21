
%% Initialization
clc;
clear
n_pair = 10000;

%% Load SyNet mat file
dsn_name = ['../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat'];
fprintf('Loading [%s] network.\n', dsn_name);
DSN_info = load(dsn_name);
n_gene = numel(DSN_info.Gene_Name);

%% Load GWAS hits over Cohort
fid = fopen('../13_GWAS_Relationship_with_DSN/DSN_iCOGS_Hits/iCOGS_Hits_Genes_MD10.0k.tsv', 'r');
% Id	-Log10(pval)	#Hit	#Hit/Size
GWAS_Info = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
fclose(fid);
if ~isequal(GWAS_Info{1}, DSN_info.Gene_Name), error(); end

%% Output makeHeatFile.py scores --heat_file
fid = fopen('./HotNet2_makeNetwork_Files/SyNet_heatfile.tsv', 'w');
for gi=1:n_gene
    fprintf(fid, '%s\t%0.2f\n', GWAS_Info{1}{gi}, GWAS_Info{2}(gi));
end
fclose(fid);

%% Output edge list: makeNetworkFiles.py --edgelist_file
fid = fopen('./HotNet2_makeNetwork_Files/SyNet_edgelist.tsv', 'w');
for pi=1:n_pair
    fprintf(fid, '%d\t%d\n', DSN_info.PP_Info(pi,1:2));
end
fclose(fid);

%% Output Index to Gene: makeNetworkFiles.py --gene_index_file
fid = fopen('./HotNet2_makeNetwork_Files/SyNet_geneindex.tsv', 'w');
for gi=1:n_gene
    fprintf(fid, '%d\t%s\n', gi, DSN_info.Gene_Name{gi});
end
fclose(fid);
