function Cur_lst = getNetNeighborsBreadthFirst(Neig_cell, Cur_lst, MAX_N_NEIGHBOR, Cur_Depth)
if numel(Cur_lst)>=MAX_N_NEIGHBOR
	Cur_lst = Cur_lst(1:MAX_N_NEIGHBOR);
elseif Cur_Depth<numel(Cur_lst)
	Cur_Depth = Cur_Depth + 1;
	Cur_lst = unique([Cur_lst Neig_cell{Cur_lst(Cur_Depth)}], 'Stable');
	Cur_lst = getNetNeighborsBreadthFirst(Neig_cell, Cur_lst, MAX_N_NEIGHBOR, Cur_Depth);
end
end