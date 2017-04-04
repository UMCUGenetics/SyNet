%% Download data
% Download RNASeq data from: https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap/HiSeqV2_percentile&host=https://tcga.xenahubs.net
% Download MicroArray data from: https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap/AgilentG4502A_07_3&host=https://tcga.xenahubs.net
% Download Clinical: https://xenabrowser.net/datapages/?dataset=TCGA.BRCA.sampleMap%2FBRCA_clinicalMatrix&host=https%3A%2F%2Ftcga.xenahubs.net
% Description of clinical data: https://gdc.cancer.gov/about-data/data-harmonization-and-generation/clinical-data-harmonization

%% Initialization
clc;
clear;
fn_expr = './data/AgilentG4502A_07_3';
fn_clic = './data/BRCA_clinicalMatrix';

%% Loading expression data
fprintf('Reading expression dataset from [%s].\n', fn_expr);
fid = fopen(fn_expr, 'r');
Patient_ID = regexp(fgetl(fid), '\t', 'split')';
Patient_ID(1) = [];
n_pat = numel(Patient_ID);
frmt_str = ['%s' repmat('%f', 1, n_pat)];
f_cell = textscan(fid, frmt_str, 'Delimiter', '\t', 'HeaderLines', 0, 'TreatAsEmpty', 'NA', 'ReturnOnError', 0, 'CollectOutput', 1);
if ~feof(fid), error(); end

Gene_Name = f_cell{1};
Gene_Expression = f_cell{2}';
[n_pat, n_gene] = size(Gene_Expression);

%% Loading clinical data
Item_Info = readtable(fn_clic, 'HeaderLines', 0, 'TreatAsEmpty', 'NA');
n_item = size(Item_Info.sampleID, 1);
Item_Map = containers.Map(Item_Info.sampleID, 1:n_item);
Patient_Info = table();
for pi=1:n_pat
	item_ind = Item_Map(Patient_ID{pi});
	Patient_Info(pi, :) = Item_Info(item_ind, :);
end
if ~isequal(Patient_ID, Patient_Info.sampleID)
	error();
end
Patient_Info.Platform = repmat({'AgilentG4502A_07_3'}, size(Patient_Info,1), 1);
Patient_Info.Subtype = repmat({'NA'}, size(Patient_Info,1), 1);

%% Replacing nans with median
for gi=1:n_gene
	is_nan = isnan(Gene_Expression(:, gi));
	if any(is_nan)
		Gene_Expression(is_nan, gi) = median(Gene_Expression(~is_nan, gi));
	end
end
if any(isnan(Gene_Expression(:))), error(); end

%% Saving Data
sav_name = 'TCGA_Combined.mat';
fprintf('Saving data in [%s]\n', sav_name);
save(sav_name, 'Gene_Expression', 'Patient_Info', 'Gene_Name');
