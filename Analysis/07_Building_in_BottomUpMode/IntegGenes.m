function [cTr, cTe] = IntegGenes(Nei_grp, zTr, lTr, zTe, Study_Index)

%% Prepare data
nTr = size(zTr, 1);
nTe = size(zTe, 1);
n_grp = numel(Nei_grp);
[~, ~, Fold_Index] = unique(Study_Index, 'Stable');
n_Fold = max(Fold_Index);

%% Combination
cTr = zeros(nTr, n_grp);
cTe = zeros(nTe, n_grp);
for gi=1:n_grp
	gene_set = Nei_grp{gi};
	n_gene = numel(gene_set);
	
	B_fold = zeros(n_gene, n_Fold);
	for fi=1:n_Fold
		B_fold(:, fi) = regress(lTr(Fold_Index==fi), zTr(Fold_Index==fi, gene_set));
	end
	B_opt = median(B_fold, 2);
	cTr(:, gi) = zTr(:, gene_set) * B_opt;
	cTe(:, gi) = zTe(:, gene_set) * B_opt;
end
end