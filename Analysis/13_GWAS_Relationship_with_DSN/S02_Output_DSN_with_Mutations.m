
%% Initialization
clc;
clear;
addpath('../../../../Useful_Sample_Codes/IsOverlapped/');
addpath('../../../../Useful_Sample_Codes/ShowProgress/');

%% Load Top Pairs
dsn_info = load('../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat');

%% Load Gene coordinates
gene_info = load('./Gene_Coordinates/Gene_Info_vGRCh37.75.mat');
n_gene = numel(gene_info.Gene_Name);

%% Load GWAS hits
fid = fopen('./iCOGS/iCOGS_SurvData_NullRem_Sorted_OnlyCoordAndPValue.tsv', 'r');
gwas_cell = textscan(fid, '%f%f%f', 10000, 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
fclose(fid);
gwas_hit = [gwas_cell{1} gwas_cell{2} gwas_cell{2} -log10(gwas_cell{3}) gwas_cell{3}];
clear gwas_cell

%% Find closest genes
MAX_DISTANCE = 10e3;
fprintf('Finding closest genes:\n');
GH_Map = containers.Map;
for gi=1:n_gene
    showprogress(gi, n_gene);
    if ~ismember(gene_info.Gene_Name{gi}, dsn_info.Gene_Name)
        continue;
    end
	dist = crdDist(gwas_hit(:,1:3), gene_info.Gene_Crd(gi,1:3));
	found_hit = find(dist<MAX_DISTANCE);
	for i=1:numel(found_hit)
		hit_ind = found_hit(i);
        if GH_Map.isKey(gene_info.Gene_Name{gi})
            GH_Map(gene_info.Gene_Name{gi}) = [GH_Map(gene_info.Gene_Name{gi}); gene_info.Gene_Crd(gi,:) gwas_hit(hit_ind,:) dist(hit_ind)];
        else
            GH_Map(gene_info.Gene_Name{gi}) = [gene_info.Gene_Crd(gi,:) gwas_hit(hit_ind,:) dist(hit_ind)];
        end
	end
end
% Hit_info = [GH_Map.values']; 
% Hit_info = vertcat(Hit_info{:});

%% Output top GWAS hits
hit_fname = sprintf('./DSN_iCOGS_Hits/iCOGS_Hits_Genes_MD%0.1fk.tsv', MAX_DISTANCE/1e3);
fprintf('Writing the GWAS hits in %s\n', hit_fname);
fid = fopen(hit_fname, 'w');
fprintf(fid, 'Id\t-Log10(pval)\t#Hit\t#Hit/Size\n');
for gi=1:numel(dsn_info.Gene_Name)
    gi_name = dsn_info.Gene_Name{gi};
    if GH_Map.isKey(gi_name)
        hit_info = GH_Map(gi_name);
        gene_score = max(hit_info(:,8));
        n_hit = size(hit_info,1);
        n_hit_norm = n_hit / (diff(hit_info(1,2:3))/1e3);
    else
        hit_info = [];
        gene_score = 0;
        n_hit = 0;
        n_hit_norm = 0;
    end
    fprintf(fid, '%s\t%0.1f\t%d\t%0.2f\n', gi_name, gene_score, n_hit, n_hit_norm);
end
fclose(fid);

%% Output pairs
SyNet_fname = sprintf('./DSN_iCOGS_Hits/DSN_SyNet_PairsWith_iCOGS_Pval_MD%0.1fk.tsv', MAX_DISTANCE/1e3);
fprintf('Writing the pair scores in %s\n', SyNet_fname);
fid = fopen(SyNet_fname, 'w');
fprintf(fid, 'Source\tTarget\tType\tWeight\tGWAS_Log10PVal\n');
n_Top = 10000;
for pi=1:n_Top
    Pair_Info = dsn_info.PP_Info(pi,:);
	gi_name = dsn_info.Gene_Name{Pair_Info(1)};
	gj_name = dsn_info.Gene_Name{Pair_Info(2)};
	gene_score = zeros(1,10);
	if GH_Map.isKey(gi_name)
		gene_score = [gene_score; GH_Map(gi_name)];
	end
	if GH_Map.isKey(gj_name)
		gene_score = [gene_score; GH_Map(gj_name)];
	end
	
	fprintf(fid, '%s\t%s\tUndirected\t%0.6f\t%0.3f\n', gi_name, gj_name, Pair_Info(15), max(gene_score(:,8)));
end
fclose(fid);

%% //// Functions
function dist = crdDist(tar_crd, ref_crd)
if size(tar_crd,2)~=3 || size(ref_crd,2)~=3, error(); end
dist = [tar_crd(:,2)-ref_crd(:,3) tar_crd(:,3)-ref_crd(:,2)];
dist = min(abs(dist), [], 2);
dist(isOL(tar_crd(:,2:3),ref_crd(:,2:3),0,0,0)) = 0;
dist(tar_crd(:,1)~=ref_crd(:,1)) = inf;
end

