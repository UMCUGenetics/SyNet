
%% Initialization
clc;
clear;
addpath('../../../../Useful_Sample_Codes/IsOverlapped/');
addpath('../../../../Useful_Sample_Codes/ShowProgress/');

%% Load Top Pairs
dsn_info = load('./Top_Pairs/Top_SyNet_AvgSynACr.mat');

%% Load Gene coordinates
gene_info = load('./Gene_Coordinates/Gene_Info_vGRCh37.75.mat');
n_gene = numel(gene_info.Gene_Name);

%% Load GWAS hits
fid = fopen('./iCOGS/iCOGS_SurvData_Sort_RemNaN.tsv', 'r');
gwas_cell = textscan(fid, '%f%f%*s%*s%*s%*s%*s%*s%*s%*s%f', 10000, 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
fclose(fid);

gwas_hit = [gwas_cell{1} gwas_cell{2} gwas_cell{2} -log10(gwas_cell{3})];
clear gwas_cell

%% Find closest genes
fprintf('Finding closest genes:\n');
GH_Map = containers.Map;
for gi=1:n_gene
    showprogress(gi, n_gene);
    if ~ismember(gene_info.Gene_Name{gi}, dsn_info.Gene_Name)
        continue;
    end
	dist = crdDist(gwas_hit(:,1:3), gene_info.Gene_Crd(gi,1:3));
	found_hit = find(dist<10e3);
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
fprintf('Writing the GWAS hits:\n');
fid = fopen('./Network_tsv/DSN_iCOGS_Hits.tsv', 'w');
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
fprintf('Writing the pair scores:\n');
fid = fopen('./Network_tsv/DSN_SyNet_vs_iCOGS.tsv', 'w');
fprintf(fid, 'Source\tTarget\tType\tWeight\n');
for pi=1:10000
	gi_name = dsn_info.Gene_Name{dsn_info.PP_Info(pi,1)};
	gj_name = dsn_info.Gene_Name{dsn_info.PP_Info(pi,2)};
	gene_score = zeros(1,9);
	if GH_Map.isKey(gi_name)
		gene_score = [gene_score; GH_Map(gi_name)];
	end
	if GH_Map.isKey(gj_name)
		gene_score = [gene_score; GH_Map(gj_name)];
	end
	
	fprintf(fid, '%s\t%s\tUndirected\t%d\n', gi_name, gj_name, max(floor(gene_score(:,8))));
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

