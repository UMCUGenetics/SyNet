function Neig_cell = getNeighborsFromAdj(Net_Adj)
if min(Net_Adj(:))<0
	error('Not implemented for negative values.');
end
n_gene = size(Net_Adj,1);
Net_Adj(1:n_gene+1:end) = 0;
[iv_mat, ii_mat] = sort(Net_Adj, 2, 'Descend');
Neig_cell = cell(n_gene, 1);
for ni=1:n_gene
	Neig_cell{ni} = [ni ii_mat(ni, iv_mat(ni,:)>0)];
end