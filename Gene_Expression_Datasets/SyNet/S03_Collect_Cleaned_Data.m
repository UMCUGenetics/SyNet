clc;
clear;

%% Use Combat (in R) to remove batch effects
% Remove_BatchEffect.R

%% Load raw data
fprintf('Loading raw data ...\n');
data_norm = load('SyNet_Normalized.mat');

%% Load batch removed data
fprintf('Loading batch removed cvs ...\n');
fid = fopen('SyNet_Normalized_Expression_Corrected.csv', 'r');
Patient_lst = regexp(fgetl(fid), '\t', 'split');
n_pat = numel(Patient_lst);
frmt_str = ['%s' repmat('%f', 1, n_pat)];
Crr_Expr = [];
Gene_Entrez = {};
batch_size = 1000;
step = 0;
while ~feof(fid)	
	fprintf('Reading next [%d] genes, step [%d] ...\n', batch_size, step);
	f_cell = textscan(fid, frmt_str, batch_size, 'Delimiter', '\t', 'HeaderLines', 0, 'ReturnOnError', 0, 'CollectOutput', 1);
	Gene_Entrez = [Gene_Entrez; f_cell{1}];
	Crr_Expr = [Crr_Expr f_cell{2}'];
	step = step + numel(f_cell{1});
end
fclose(fid);
Gene_Entrez = strrep(Gene_Entrez, '"', '');
if ~isequal(Gene_Entrez, data_norm.Gene_Entrez), error(); end

%% Converting Entrez to Hugo GeneName
% ! wget ftp://ftp.ncbi.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz
% # For complete description refer to ftp://ftp.ncbi.nih.gov/gene/DATA/README
% ! zcat Homo_sapiens.gene_info | awk '{print($2"\t"$3"\t"$5)}' > Entrez2GeneName.txt
fid = fopen('./Entrez2GeneName.txt', 'r');
f_cell = textscan(fid, '%s%s%*s', 'Delimiter', '\t', 'HeaderLines', 1, 'ReturnOnError', 0);
if ~feof(fid), error(); end
n_gene = numel(Gene_Entrez);
Gene_Name = cell(n_gene, 1);
Ent2Name = containers.Map(f_cell{1}, f_cell{2});
for gi=1:n_gene
	if Ent2Name.isKey(Gene_Entrez{gi})
		Gene_Name{gi} = Ent2Name(Gene_Entrez{gi});
	end
end

%% Indexify studies
[Study_Name, ~, Study_Index] = unique(data_norm.Patient_Info.StudyName, 'Stable');
Study_Name = strcat(Study_Name, [repmat({'-ACES'},12,1); repmat({'-HAIBE'},36,1); {''}; {''}]);

%% Saving data
data_crr.Gene_Expression = Crr_Expr;
data_crr.Gene_Entrez = Gene_Entrez;
data_crr.Patient_Info = data_norm.Patient_Info;
data_crr.Patient_Label = double(data_norm.Patient_Info.Prognostic_Status);
data_crr.Gene_Name = Gene_Name;
data_crr.Study_Name = Study_Name;
data_crr.Study_Index = Study_Index;
sav_name = 'SyNet_BatchCorrected.mat';
fprintf('Saving corrected data in [%s]\n', sav_name);
save(sav_name, '-struct', 'data_crr');
