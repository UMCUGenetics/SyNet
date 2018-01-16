clc;
clear

%% Initialization
addpath('../17_Frequency_of_SyNetLinks_in_Biological_Netwrok/');
addpath('../_Utilities/');
MAX_DISTANCE = 10000;
MAX_SyNet_Pairs = 03544;
%MAX_SyNet_Pairs = 50000;
% MAX_N_SNP = 504752;
% MAX_N_SNP = 53837; % pval = 0.00499986
% MAX_N_SNP = 10962; % pval = 0.000999939
% MAX_N_SNP = 10000; % pval = 0.000999939
MAX_N_SNP = 1104; % pval = 0.0000999498
% MAX_N_SNP = 1000; % pval = 0.0000910276
% MAX_N_SNP = 500; % pval = 0.0000501868
Ref_Name = 'AvgSynACr';

%% Load SyNet mat file
dsn_name = ['../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet_' Ref_Name '.mat'];
fprintf('Loading [%s] network.\n', dsn_name);
DSN_info = load(dsn_name);
n_gene = numel(DSN_info.Gene_Name);

%% Load GWAS hits over Cohort
gwas_name = sprintf('../13_GWAS_Relationship_with_DSN/iCOGS_Hits/iCOGS_Hits_Genes_NTS%0.1fk_MD%0.1fk.tsv', MAX_N_SNP/1e3, MAX_DISTANCE/1e3);
fprintf('Loading GWAS data from [%s]\n', gwas_name);
fid = fopen(gwas_name, 'r');
% Id	-Log10(pval)	#Hit	#Hit/Size
GWAS_Info = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
fclose(fid);
if ~isequal(GWAS_Info{1}, DSN_info.Gene_Name), error(); end

%% Main loop
for ri=1:2
    
    %% Source preparation
    if ri==1
        Pair_GName = DSN_info.Gene_Name;
    else
        rind = randperm(n_gene);
        Pair_GName = DSN_info.Gene_Name(rind);
        Ref_Name = [Ref_Name '-SF'];
    end
    Src_dir = sprintf('./HotNet2_Files/RF-%s_MNP-%06d_MNS-%06d_MD-%0.0fk/', Ref_Name, MAX_SyNet_Pairs, MAX_N_SNP, MAX_DISTANCE/1e3);
    if exist(Src_dir, 'dir')
        rmdir(Src_dir, 's');
    end
    fprintf('Saving input files in [%s]\n', Src_dir);
    mkdir([Src_dir 'TMP_DIR']);
    mkdir([Src_dir 'Input_DIR']);
    mkdir([Src_dir 'Output_DIR']);
    
    %% Output makeHeatFile.py pvalues --heat_file
    GWAS_Info{2} = GWAS_Info{2} / max(GWAS_Info{2});
    fid = fopen([Src_dir 'Input_DIR/iCOGS_LogPval_Normalized.tsv'], 'w');
    for gi=1:n_gene
        fprintf(fid, '%s\t%0.5f\n', GWAS_Info{1}{gi}, GWAS_Info{2}(gi));
    end
    fclose(fid);
    
    %% Output makeHeatFile.py # Hits notmalized by size --heat_file
    GWAS_Info{4} = GWAS_Info{4} / max(GWAS_Info{4});
    fid = fopen([Src_dir 'Input_DIR/iCOGS_NHitSize_Normalized.tsv'], 'w');
    for gi=1:n_gene
        fprintf(fid, '%s\t%0.5f\n', GWAS_Info{1}{gi}, GWAS_Info{4}(gi));
    end
    fclose(fid);
    
    %% Output edge list: makeNetworkFiles.py --edgelist_file
    Pair_List = DSN_info.PP_Info(1:MAX_SyNet_Pairs,1:2);
    fid = fopen([Src_dir 'Input_DIR/EdgeList.tsv'], 'w');
    for pi=1:MAX_SyNet_Pairs
        fprintf(fid, '%d\t%d\n', Pair_List(pi,1:2));
    end
    fclose(fid);
    
    %% Output Index to Gene: makeNetworkFiles.py --gene_index_file
    fid = fopen([Src_dir 'Input_DIR/GeneIndex.tsv'], 'w');
    for gi=1:n_gene
        fprintf(fid, '%d\t%s\n', gi, Pair_GName{gi});
    end
    fclose(fid);
end

