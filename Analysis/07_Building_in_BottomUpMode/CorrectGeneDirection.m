function [xTr, xTe] = CorrectGeneDirection(xTr, xTe, lTr)
n_gene = size(xTr, 2);
gene_dir = corr(xTr, lTr, 'Type', 'Spearman');
for gi=1:n_gene
	if gene_dir(gi)<0
		xTr(:, gi) = -xTr(:, gi);
		xTe(:, gi) = -xTe(:, gi);
	end
end
end
