function [patch_h, bin_freq] = ViolinEx(x_Pos, edge_lst, dist, opt)
if ~exist('opt', 'var'), opt = struct; end
if ~isfield(opt, 'BarColor')
    opt.BarColor = [0.1 0.1 1];
end
n_bin = numel(edge_lst)-1;
hold on

%% Make PDF
bin_freq = zeros(1, n_bin);
for bi=1:n_bin
    if bi==n_bin
        has_ol = edge_lst(bi)<=dist & dist <= edge_lst(bi+1);
    else
        has_ol = edge_lst(bi)<=dist & dist <  edge_lst(bi+1);
    end
    bin_freq(bi) = sum(has_ol);
end

%% Normalize and adjustments
if ~isfield(opt, 'Normalize') || opt.Normalize==1
    bin_nrm = bin_freq/max(bin_freq);
else
    bin_nrm = bin_freq;
end
if isfield(opt, 'ShrinkFactor')
    bin_nrm = bin_nrm * opt.ShrinkFactor;
end
if isfield(opt, 'Reverse') && opt.Reverse==1
    bin_nrm = - bin_nrm;
end

%% Plot bars
for bi=1:n_bin
    patch_h(bi,1) = patch([x_Pos x_Pos+bin_nrm([bi bi]) x_Pos], edge_lst([bi bi bi+1 bi+1]), 'b', 'EdgeColor', 'None', 'FaceColor', opt.BarColor);
end
end
