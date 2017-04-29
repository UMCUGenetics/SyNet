function pred = nmc(x_tr, l_tr, x_te)
Class_list = unique(l_tr);
n_class = numel(Class_list);
mean_tr = zeros(n_class, size(x_tr, 2));
for ci=1:n_class
    mean_tr(ci,:) = mean(x_tr(l_tr==Class_list(ci), :), 1);
end
di = pdist2(mean_tr, x_te);
[~, min_ind] = min(di, [], 1);
pred = Class_list(min_ind);
end