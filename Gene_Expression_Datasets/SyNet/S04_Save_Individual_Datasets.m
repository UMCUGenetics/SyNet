clc;
clear;

%% Load data
fprintf('Loading batch corrected data ...\n');
data_synet = load('SyNet_BatchCorrected.mat');

%% Save data per source
source_lst = {'ACES', 'METABRIC', 'TCGA'};
for si=1:numel(source_lst)
	src_name = source_lst{si};
	sav_name = ['../' src_name '/' src_name '_BatchCorrected.mat'];
	fprintf('Saving [%s] in [%s] ...\n', src_name, sav_name);
	
	is_in = ismember(data_synet.Patient_Info.Source_Study, src_name);
	out.Gene_Entrez = data_synet.Gene_Entrez;
	out.Gene_Name = data_synet.Gene_Name;
	out.Gene_Expression = data_synet.Gene_Expression(is_in, :);
	out.Patient_Info = data_synet.Patient_Info(is_in, :);
	out.Patient_Label = data_synet.Patient_Label(is_in);
	out.Study_Index = data_synet.Study_Index(is_in);
	[out.Study_Name, ~, out.Study_Index] = unique(data_synet.Study_Name(data_synet.Study_Index(is_in)), 'stable');
	
	save(sav_name, '-struct', 'out');
	clear out
end
