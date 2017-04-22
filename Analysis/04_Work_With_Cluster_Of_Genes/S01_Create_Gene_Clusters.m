%% Initialization
clc;
clear
addpath('../../../../Useful_Sample_Codes/ShowProgress/');

%% Load data
GeneExpression_Path = '../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat';
load(GeneExpression_Path, 'Gene_Expression', 'Patient_Label', 'Gene_Name', 'Study_Index', 'Patient_Info');
zData = quantilenorm(Gene_Expression);
n_gene = size(Gene_Expression, 2);
SurvivalTime = Patient_Info.SurvivalTime;

%% Plot cluster size
% Z = linkage(zData', 'average', 'correlation');
% Cls_ind = cluster(Z, 'cutoff', 0.6);
% freq = accumarray(Cls_ind, 1);
% plot(sort(freq));

%% Calculate correlation
cr_mat = corr(Gene_Expression, 'Type', 'Spearman');

%% Sorting of neighbors
[cr_val, cr_ind] = sort(abs(cr_mat), 2, 'Descend');

%% Grouping
nei_grp = cell(n_gene, 1);
is_used = false(n_gene, 1);
for gi=1:n_gene
	showprogress(gi, n_gene);
	nei_grp{gi} = cr_ind(gi, cr_val(gi,:)>=0.5);
	if numel(nei_grp{gi})>1
		is_used(nei_grp{gi}) = 1;
	end
end

%% Measure size
grp_size = cellfun('length', nei_grp);
[~, sid] = sort(grp_size, 'Descend');
nei_grp = nei_grp(sid);
grp_size = cellfun('length', nei_grp);
plot(grp_size);

%% Identify clustered genes
has_grp = unique([nei_grp{grp_size>1}]);
cr_grp = corr(Gene_Expression(:, has_grp), SurvivalTime, 'Type', 'Spearman', 'rows', 'complete');

no_grp = setdiff(1:n_gene, has_grp);
cr_ngp = corr(Gene_Expression(:, no_grp), SurvivalTime, 'Type', 'Spearman', 'rows', 'complete');

close all
hold on
plot(sort(abs(cr_grp), 'Descend'), 'b', 'LineWidth', 2);
% plot(sort(abs(cr_grp(randperm(numel(cr_grp), numel(cr_ngp)))), 'Descend'), 'b');
% plot(sort(abs(cr_ngp), 'Descend'), 'r');
plot(sort(abs(cr_ngp(randperm(numel(cr_ngp), numel(cr_grp)))), 'Descend'), 'r', 'LineWidth', 2);
legend({'Grouped' 'Individuals'}, 'FontWeight', 'Bold', 'FontSize', 14);

%% Count how many are relevant
load('../01_Pairwise_Evaluation_of_Genes/Top_Pairs_CV01/Top_SyNet.mat', 'PP_Info');
tresh_lst = [1e2 5e2 1e3 5e3 10e3 20e3 50e3 100e3];
n_tresh = numel(tresh_lst);
top_freq = zeros(n_tresh, 5);
for ti=1:n_tresh
	Top_Pairs = unique(PP_Info(1:tresh_lst(ti), 1:2));
	
	top_freq(ti,1) = sum(ismember(has_grp, Top_Pairs));
	top_rnd = randperm(n_gene, numel(has_grp));
	top_freq(ti,2) = sum(ismember( top_rnd, Top_Pairs));
	
	top_freq(ti,3) = sum(ismember( no_grp, Top_Pairs));
	top_rnd = randperm(n_gene, numel(no_grp));
	top_freq(ti,4) = sum(ismember( top_rnd, Top_Pairs));
	
	top_freq(ti,5) = numel(Top_Pairs);
end

figure('Position', [100 100 1400 500]);
bar(top_freq(:,1:4)*100./repmat(top_freq(:,5),1,4), 'Grouped');
ylim([0 100]);
ylabel('Frequency (%)');
xlabel('#top pairs selected from SPADE');
set(gca, 'XTick', 1:n_tresh, 'XTickLabel', tresh_lst, 'FontWeight', 'Bold', 'FontSize', 14);

colormap([
	0.2 0.2 1.0
	0.6 0.6 0.9
	1.0 0.2 0.2
	0.9 0.6 0.6
	]);
legend({'Grouped', 'Random (grp)', 'Individuals', 'Random (indv)'}, 'FontWeight', 'Bold', 'FontSize', 14, 'Location', 'NorthEastOutside');
for ti=1:n_tresh
	text(ti, 95, sprintf('#Genes=%d', top_freq(ti,5)), 'FontWeight', 'Bold', 'VerticalAlignment', 'Middle', 'HorizontalAlignment', 'Center', 'FontSize', 12);
end

%% Select top correlating genes
n_top = 100000;
[cr_val, cr_ind] = sort(abs(cr_mat(:)), 'Descend');
cr_ind = cr_ind(1:n_top);
cr_val = cr_val(1:n_top);

[Pair_Corr(:,1), Pair_Corr(:,2)] = ind2sub([n_gene, n_gene], cr_ind);
Pair_Corr(:, 3) = cr_val;

gi = 1;
close all
hold on
lbl = Patient_Label;
% lbl = Study_Index;
for l=unique(lbl(:))'
	plot(Gene_Expression(lbl==l,Pair_Corr(gi,1)), Gene_Expression(lbl==l,Pair_Corr(gi,2)), '.');
end

%% Compare with Census
fid = fopen('../../../../../Dataset/Census/cancer_gene_census_v8.0_gene_symbol.tsv');
census_gn = textscan(fid, '%s', 'Delimiter', '\t', 'ReturnOnError', 0, 'HeaderLine', 0);
if ~feof(fid), error(); end
fclose(fid);
census_gn = unique(census_gn{1});

cen_freq(1,1) = sum(ismember(Gene_Name(has_grp), census_gn));
top_rnd = randperm(n_gene, numel(has_grp));
cen_freq(1,2) = sum(ismember(Gene_Name(top_rnd), census_gn));
	
cen_freq(1,3) = sum(ismember(Gene_Name(no_grp), census_gn));
top_rnd = randperm(n_gene, numel(no_grp));
cen_freq(1,4) = sum(ismember(Gene_Name(top_rnd), census_gn));

cen_freq(1,5) = numel(census_gn);

figure('Position', [100 100 1400 500]);
bar(cen_freq(1:4)*100./cen_freq(:,5), 'Grouped');
ylim([0 100]);
ylabel('Frequency (%)');
xlabel('#top pairs selected from SPADE');
set(gca, 'XTick', 1:n_tresh, 'XTickLabel', {'Grouped', 'Random (grp)', 'Individuals', 'Random (indv)'}, 'FontWeight', 'Bold', 'FontSize', 14);

%% Compute correlation of top pairs
load('/Users/amin/Technical/My_Works/Deep_Learning/113_Organize_SPADE_Codes/Analysis/01_Pairwise_Evaluation_of_Genes/Top_Pairs_CV01/Top_SyNet.mat')
n_top = 5000;
tp = PP_Info(1:n_top,:);
tp_cr = zeros(n_top,1);
for pi=1:n_top
	tp_cr(pi,1) = corr(zData(:, tp(pi,1)), zData(:, tp(pi,2)), 'Type', 'Spearman');
end
cr_nrm = floor(abs(tp_cr)*7)+1;
scatter(tp(:,6), tp(:,7), cr_nrm, cr_nrm);
colormap(jet(7));

