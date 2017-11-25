function Ind_List = UnifyGeneList(ref_lst, que_lst)
n_ref = numel(ref_lst);
n_que = numel(que_lst);
gMap = containers.Map(que_lst, 1:n_que);

Ind_List = zeros(n_ref, 1);
for gi=1:n_ref
    Ind_List(gi) = gMap(ref_lst{gi});
end
if ~isequal(ref_lst, que_lst(Ind_List)), error(); end
end

