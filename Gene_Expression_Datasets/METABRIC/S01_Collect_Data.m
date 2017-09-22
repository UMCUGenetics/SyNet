%% Download data
% Download data from: https://www.synapse.org/#!Synapse:syn2133309

%% Initialization
clc;
clear;
fn_expr = './data/Metabric_AllProbs.mat';
fn_clic = './data/Metabric_Clinical.csv';

%% Loading expression data
fprintf('Reading expression dataset from [%s].\n', fn_expr);
data_expr = load(fn_expr);
All_Expression = data_expr.ExpressionData;
[n_pat, n_prob] = size(All_Expression);
Pat_Map = containers.Map(data_expr.SampleID, 1:n_pat);
Prob_ID = data_expr.Entrez_ID;

%% Loading clinical data
Patient_Info = readtable(fn_clic, 'HeaderLines', 0, 'TreatAsEmpty', 'NA');
n_clic = size(Patient_Info,1);
if ~isequal(data_expr.SampleID, {Patient_Info.patient_id{:}}')
	error();
end
Patient_Info.Platform = repmat({'Illumina_Human_WG-v3'}, size(Patient_Info,1), 1);

%% Combine probs to Entrez
[Gene_Entrez, ~, Prob_idx] = unique(Prob_ID);
Prob_grp = accumarray(Prob_idx, 1:numel(Prob_idx), [], @(i) {i});
n_gene = numel(Gene_Entrez);
Gene_Expression = nan(n_pat, n_gene);
for gi=1:n_gene
	prb_set = Prob_grp{gi};
	if ~isequal(Prob_ID{prb_set(1)}, Gene_Entrez{gi}), error(); end
	
	prb_std = std(All_Expression(:, prb_set));
	if any(isnan(prb_std)), error('Std of a prob is NaN.'); end
	[~, sid] = sort(prb_std, 'Descend');
	prb_set = prb_set(sid);
	Gene_Expression(:, gi) = All_Expression(:, prb_set(1));
end
if any(isnan(Gene_Expression(:))), error(); end

%% Saving Data
sav_name = 'METABRIC_Combined.mat';
fprintf('Saving data in [%s]\n', sav_name);
save(sav_name, 'Gene_Expression', 'Patient_Info', 'Gene_Entrez');
