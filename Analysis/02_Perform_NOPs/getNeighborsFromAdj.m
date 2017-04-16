function Neig_cell = getNeighborsFromAdj(Net_Adj)
if min(Net_Adj(:))<0
	error('Not implemented for negative values.');
end
n_gene = size(Net_Adj,1);

%% Add eps to binary networks to avoid selecting the same neighbors according to their indices and not based on interaction strength
if isequal(unique(Net_Adj(:)), [0; 1])
	fprintf('[i] Adjacency matrix is binary, adding eps to non-zero elements.\n');
	is_zero = Net_Adj(:) == 0;
	Net_Adj = Net_Adj + rand(n_gene)*eps;
	Net_Adj(is_zero) = 0;
end

%% Neighborhood selection
Net_Adj(1:n_gene+1:end) = 0;
[iv_mat, ii_mat] = sort(Net_Adj, 2, 'Descend');
Neig_cell = cell(n_gene, 1);
for ni=1:n_gene
	Neig_cell{ni} = [ni ii_mat(ni, iv_mat(ni,:)>0)];
end
