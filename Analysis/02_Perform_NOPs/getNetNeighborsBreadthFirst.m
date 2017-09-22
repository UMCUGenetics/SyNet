function Cur_lst = getNetNeighborsBreadthFirst(Neig_cell, Cur_lst, n_neighbor, Depth_ind)
if numel(Cur_lst)>=n_neighbor
	Cur_lst = Cur_lst(1:n_neighbor);
elseif Depth_ind<numel(Cur_lst)
	Depth_ind = Depth_ind + 1;
	Cur_lst = unique([Cur_lst Neig_cell{Cur_lst(Depth_ind)}], 'stable');
	Cur_lst = getNetNeighborsBreadthFirst(Neig_cell, Cur_lst, n_neighbor, Depth_ind);
end
end