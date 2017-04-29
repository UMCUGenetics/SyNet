function [cTr, cTe] = CombineGenes(Nei_grp, zTr, zTe, MAX_N_NEIGHBORS)

%% Prepare data
nTr = size(zTr, 1);
nTe = size(zTe, 1);
n_grp = numel(Nei_grp);

%% Combination
cTr = zeros(nTr, n_grp);
cTe = zeros(nTe, n_grp);
for gi=1:n_grp
	gene_set = Nei_grp{gi};
	if numel(gene_set)>MAX_N_NEIGHBORS
		gene_set = gene_set(1:MAX_N_NEIGHBORS);
	end
	
	cTr(:, gi) = mean(zTr(:, gene_set), 2);
	cTe(:, gi) = mean(zTe(:, gene_set), 2);
end
end