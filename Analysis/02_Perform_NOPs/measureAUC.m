function res = measureAUC(x, l, n_rep)
if ~isequal(size(x), size(l)), error(); end

n_pop = numel(l);
n_sample = floor(n_pop*0.7);
auc_lst = zeros(n_rep, 1);
for ri=1:n_rep
	ind = randperm(n_pop, n_sample);
	auc_lst(ri) = getAUC(l(ind), x(ind), 50);
end
res = median(auc_lst);
end