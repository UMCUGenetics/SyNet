clc;
clear;

%% Load GWAS hits over Cohort
fid = fopen('./DSN_iCOGS_Hits/iCOGS_Hits_Genes_MD10.0k.tsv', 'r');
% Id	-Log10(pval)	#Hit	#Hit/Size
gwas_hits = textscan(fid, '%s%f%f%f', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
fclose(fid);
Gene_Name = gwas_hits{1};
gwas_hits = [gwas_hits{2} gwas_hits{3} gwas_hits{4}];
n_gene = numel(Gene_Name);

%% Load SyNet pairs
n_top = 10000;
Score_ind = 15;
dsn_info = load('../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat');
if ~isequal(Gene_Name, dsn_info.Gene_Name), error(); end
% All_Info = [dsn_info.PP_Info(1:n_top, [1 2 Score_ind]); dsn_info.NP_Info(1:n_top, [1 2 Score_ind])];
All_Info = dsn_info.PP_Info(1:n_top, [1 2 15]);

%% Assign pair fitness to gene fitness
IG_Info = zeros(n_gene, 1);
for pi=size(All_Info,1):-1:1
    IG_Info(All_Info(pi,1:2)) = All_Info(pi,3);
end
gwas_hits(:,4) = IG_Info;

%% Categorize and visualize genes
is_Hit = gwas_hits(:,2)>0;
[h, pval, ci, stats] = ttest2(IG_Info(is_Hit), IG_Info(~is_Hit));

close all
histogram(IG_Info(~is_Hit), 'FaceColor', [0.8 0.8 0.8], 'Normalization', 'pdf', 'BinLimits', [0 1], 'NumBins', 20);
hold on
histogram(IG_Info(is_Hit), 'FaceColor', 'g', 'Normalization', 'pdf', 'BinLimits', [0 1], 'NumBins', 20);
legend({'Others' 'GWas Hit'});
title(sprintf('#Pairs=%d, #Hits=%d, p-value=%0.4f', size(All_Info,1), sum(is_Hit), pval), 'FontSize', 12);
corr(gwas_hits, 'type', 'Spearman')


