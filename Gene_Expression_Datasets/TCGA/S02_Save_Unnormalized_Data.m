clc;
clear;

%% Load un-normalized data
synet_data = load('../../Gene_Expression_Datasets/SyNet/SyNet_Normalized_Par.mat');

%% Select TCGA
is_in = ismember(synet_data.Patient_Info.Source_Study, 'TCGA') & ~isnan(synet_data.Patient_Info.Prognostic_Status);
Gene_Expression = synet_data.Gene_Expression(is_in,:);
Patient_Info = synet_data.Patient_Info(is_in, :);
Patient_Label = double(Patient_Info.Prognostic_Status);
Gene_Entrez = synet_data.Gene_Entrez;
[Study_Name, ~, Study_Index] = unique(Patient_Info.StudyName, 'Stable');

%% Converting Entrez to Hugo GeneName
fid = fopen('../SyNet/Entrez2GeneName.txt', 'r');
f_cell = textscan(fid, '%s%s%*s', 'Delimiter', '\t', 'HeaderLines', 1, 'ReturnOnError', 0);
if ~feof(fid), error(); end
fclose(fid);
n_gene = numel(Gene_Entrez);
Gene_Name = cell(n_gene, 1);
Ent2Name = containers.Map(f_cell{1}, f_cell{2});
for gi=1:n_gene
	if Ent2Name.isKey(Gene_Entrez{gi})
		Gene_Name{gi} = Ent2Name(Gene_Entrez{gi});
	else
		error();
	end
end

%% Save data
save('./TCGA_Unnormalized.mat', 'Gene_Expression', 'Patient_Info', 'Patient_Label', ...
	'Gene_Name', 'Gene_Entrez', 'Study_Name', 'Study_Index');