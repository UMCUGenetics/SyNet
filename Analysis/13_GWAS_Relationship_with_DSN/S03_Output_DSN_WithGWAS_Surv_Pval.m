
%% Initialization
clear;
addpath('../../../../Useful_Sample_Codes/IsOverlapped/');

%% Load Gene coordinates
gene_info = load('./Gene_Coordinates/Gene_Info_vGRCh37.75.mat');
n_gene = numel(gene_info.Gene_Name);

%% Load GWAS hit
fid = fopen('./GWAS_data/79_breast_cancer_susceptibility_loci_previously_reported_in_studies.tsv', 'r');
gwasOld_cell = textscan(fid, '%*s%f%f%*s%*s%*s%*s%*s%*s%*s%*s%f', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
if ~feof(fid), error(); end
fclose(fid);

fid = fopen('./GWAS_data/15_breast_cancer_susceptibility_loci_new.tsv', 'r');
gwasNew_cell = textscan(fid, '%*s%f%f%*s%*s%*s%*s%*s%*s', 'HeaderLines', 1, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
if ~feof(fid), error(); end
fclose(fid);

hit_scr = -log10(gwasOld_cell{3});
gwas_hit = [
	gwasOld_cell{1} gwasOld_cell{2} gwasOld_cell{2} hit_scr
	gwasNew_cell{1} gwasNew_cell{2} gwasNew_cell{2} repmat(max(hit_scr)+1, numel(gwasNew_cell{2}), 1)
];

%% Sorting GWAS hits
[~, sid] = sort(gwas_hit(:,4), 'Descend');
gwas_hit = gwas_hit(sid,:);
n_hit = size(gwas_hit, 1);

%% Find closest genes
GH_Map = containers.Map;
for hi=1:n_hit
	dist = crdDist(gene_info.Gene_Crd(:,1:3), gwas_hit(hi,1:3));
	has_ovl = find(dist<50e3);
	for i=1:numel(has_ovl)
		gi = has_ovl(i);
		GH_Map(gene_info.Gene_Name{gi}) = [gene_info.Gene_Crd(gi,:) gwas_hit(hi,:) dist(gi)];
	end
end
Hit_info = [GH_Map.values']; 
Hit_info = vertcat(Hit_info{:});

%% Load Top Pairs
top_info = load('./Top_Pairs/Top_SyNet_AvgSynACr.mat');

%% Output pairs
fid = fopen('./Network_tsv/DSN_SyNet_vs_GWAS.tsv', 'w');
for pi=1:10000
	gi_name = top_info.Gene_Name{top_info.PP_Info(pi,1)};
	gj_name = top_info.Gene_Name{top_info.PP_Info(pi,2)};
	gene_score = zeros(2,9);
	if GH_Map.isKey(gi_name)
		gene_score(1,:) = GH_Map(gi_name);
	end
	if GH_Map.isKey(gj_name)
		gene_score(2,:) = GH_Map(gj_name);
	end
	
	fprintf(fid, '%s\t%s\t%d\n', gi_name, gj_name, max(floor(gene_score(:,8))));
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

