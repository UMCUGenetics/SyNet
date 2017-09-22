clc;
clear;

%% Load metabric
meta_data = load('METABRIC_Combined.mat');

%% Converting Entrez to Hugo GeneName
% ! wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
% # For complete description refer to ftp://ftp.ncbi.nih.gov/gene/DATA/README
% ! zcat Homo_sapiens.gene_info | awk '{print($2"\t"$3"\t"$5)}' > Entrez2GeneName.txt
fid = fopen('../SyNet/Entrez2GeneName.txt', 'r');
f_cell = textscan(fid, '%s%s%*s', 'Delimiter', '\t', 'HeaderLines', 1, 'ReturnOnError', 0);
if ~feof(fid), error(); end
n_gene = numel(meta_data.Gene_Entrez);
meta_data.Gene_Name = cell(n_gene, 1);
Ent2Name = containers.Map(f_cell{1}, f_cell{2});
for gi=1:n_gene
	if Ent2Name.isKey(meta_data.Gene_Entrez{gi})
		meta_data.Gene_Name{gi} = Ent2Name(meta_data.Gene_Entrez{gi});
	else
		meta_data.Gene_Name{gi} = '';
	end
end

%% Saving data
save('METABRIC_Combined.mat', '-struct', 'meta_data');