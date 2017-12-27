function Neig_cell = getNeighborsFromAdj(Net_Adj, MAX_NEIGH_SIZE)
if ~exist('MAX_NEIGH_SIZE', 'var')
	MAX_NEIGH_SIZE = 200;
end
if min(Net_Adj(:))<0
	error('Not implemented for negative values.');
end
n_gene = size(Net_Adj,1);

%% Add eps to binary networks to avoid selecting the same neighbors according to their indices and not based on interaction strength
% fprintf('[i] Adding eps to non-zero elements to perturb the equal interaction strength.\n');
is_zero = Net_Adj(:) == 0;
Net_Adj = Net_Adj + rand(n_gene)*eps;
Net_Adj(is_zero) = 0;

%% Neighborhood selection
Net_Adj(1:n_gene+1:end) = 0;
[iv_mat, ii_mat] = sort(Net_Adj, 2, 'Descend');
Neig_cell = cell(n_gene, 1);
for gi=1:n_gene
	gene_nei = ii_mat(gi, iv_mat(gi,:)>0);
	if numel(gene_nei)>MAX_NEIGH_SIZE
		gene_nei = gene_nei(1:MAX_NEIGH_SIZE);
	end
	Neig_cell{gi} = [gi gene_nei];
end
end