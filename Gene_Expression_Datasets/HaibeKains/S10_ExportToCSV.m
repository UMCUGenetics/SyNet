clc;
clear;

%% Initialize
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
load('./HaibeKains_Combined.mat');
n_patient = size(GeneExpression, 1);
[Study_Name, ~, Study_Index] = unique({Patient_Info.StudyName}, 'Stable');
n_study = numel(Study_Name);

%% Remove probs with many nans
fprintf('Discarding probs/patients with many #NANs ...\n');
is_nan = isnan(GeneExpression);
is_val = sum(~is_nan,1)/n_patient > 0.4;
GeneExpression = GeneExpression(:, is_val);
Gene_Entrez = Gene_Entrez(is_val);
Gene_Name = Gene_Name(is_val);
Prob_ID = Prob_ID(is_val);
n_gene = size(GeneExpression, 2);

%% Normalize data (ignore NANs)
fprintf('Normalizing probs: ');
is_nan = isnan(GeneExpression);
for gi=1:n_gene
	showprogress(gi, n_gene);
	for si=1:n_study
		is_good = Study_Index==si & ~is_nan(:, gi);
		is_bad  = Study_Index==si &  is_nan(:, gi);
		GeneExpression(is_good, gi) = zscore(GeneExpression(is_good, gi));
		GeneExpression( is_bad, gi) = median(GeneExpression(is_good, gi));
	end
end

%% Save the data
save('./HaibeKains_NaNFiltered.mat', 'GeneExpression', 'Gene_Entrez', 'Gene_Name', 'Prob_ID', 'Patient_Info');


