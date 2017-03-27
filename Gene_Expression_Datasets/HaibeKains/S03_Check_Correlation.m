clc;
clear;

%% Initialize
load('./HaibeKains_Combined.mat', 'GeneExpression', 'Patient_Info');
[Study_Name, ~, Study_Index] = unique({Patient_Info.StudyName}, 'Stable');
n_study = max(Study_Index);
[n_patient, n_gene] = size(GeneExpression);

%% Remove genes with many nans
del_ind = [];
for gi=1:n_gene
	is_nan = isnan(GeneExpression(:, gi));
	if sum(is_nan)/n_patient>0.6
		del_ind = [del_ind gi];
	else
		GeneExpression(is_nan, gi) = median(GeneExpression(~is_nan, gi));
	end
end
GeneExpression(:, del_ind) = [];
n_gene = size(GeneExpression, 2);

%% Main loop
for si=1:n_study
	in = Study_Index == si;
	GeneExpression(in, :) = zscore(GeneExpression(in,:));
end

%% Calculate correlation
cr = corr(GeneExpression', 'Type', 'Spearman');
cr(1:n_patient+1:end) = 0;
has_cr = max(abs(cr))>0.8;
img = cr(has_cr, has_cr);

%% Plotting
n_item = size(img, 1);
imagesc(img);
colormap(jet(10));

set(gca, 'XTick', 1:n_item, 'XTickLabel', {Patient_Info(has_cr).StudyName}, 'XTickLabelRotation', 45, ...
		 'YTick', 1:n_item, 'YTickLabel', {Patient_Info(has_cr).StudyName});

