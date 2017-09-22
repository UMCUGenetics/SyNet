function [Freq_mat, Xedges, Yedges, XClass, YClass] = plotAxisCvg(Pair_lst, n_xbin, n_ybin)

%% Discretization
x_edge = linspace(min(Pair_lst(:,1)), max(Pair_lst(:,1)), n_xbin+1);
y_edge = linspace(min(Pair_lst(:,2)), max(Pair_lst(:,2)), n_ybin+1);

%% Counting
[Freq_mat, Xedges, Yedges, XClass, YClass] = histcounts2(Pair_lst(:,1), Pair_lst(:,2), x_edge, y_edge);
end