function [pred_lbl, opt_K, dist] = knn(xTr, lTr, xTe, opt_info)
%% Initialization
n_Te = size(xTe, 1);
[class_list, ~, class_uid] = unique(lTr);
n_class = max(class_uid);
if ~exist('opt_info', 'var')
    opt_info = struct();
end
if isfield(opt_info, 'WeightStratified')
    class_freq = histcounts(class_uid, 1:n_class+1);
    weights = 1./class_freq;
    weights = weights/sum(weights);
else
    weights = repmat(1/n_class, 1, n_class);
end

%% Use or find K
if isfield(opt_info, 'iCvPar')
    fprintf('Attempting to identify the best K ...\n');
    K_lst = 1:2:10;
    n_K = numel(K_lst);
    [~, ~, iCvPar] = unique(opt_info.iCvPar, 'Stable');
    in_opt = rmfield(opt_info, 'iCvPar');
    n_fold = max(iCvPar);
	in_auc = zeros(n_fold, n_K);
    for fi=1:n_fold
        fprintf('[%02d/%02d] AUC is: ', fi, n_fold);
        iTr = iCvPar~=fi;
        iTe = iCvPar==fi;
        in_opt.('TeTrDistance') = pdist2(xTr(iTe,:), xTr(iTr,:));
        for ki=1:n_K
            pred_lbl = knn(xTr(iTr,:), lTr(iTr,:), xTr(iTe,:), setfield(in_opt, 'K', K_lst(ki)));
            in_auc(fi, ki) = getAUC(lTr(iTe), pred_lbl, 50);
            fprintf(', K%d=%0.2f', K_lst(ki), in_auc(fi,ki)*100);
        end
        fprintf('\n');
    end
    fprintf('Final AUC is:     ');
    disp(mean(in_auc,1)*100);
    [~, opt_Kind] = max(mean(in_auc,1));
    opt_K = K_lst(opt_Kind);
    fprintf('Cross-validation is finished, best K is: [%d]\n', opt_K);
else
    opt_K = opt_info.K;
end

%% Compute distance
if isfield(opt_info, 'TeTrDistance')
    dist = opt_info.TeTrDistance;
else
    dist = pdist2(xTe, xTr);
end
[~, nei_ind] = sort(dist, 2, 'ascend');

%% Find top neighbors
pred_lbl = zeros(n_Te, 1);
for si=1:n_Te
    pred_set = class_uid(nei_ind(si,1:opt_K));
    pred_freq = histcounts(pred_set, 1:n_class+1);
    pred_freq = pred_freq/sum(pred_freq) .* weights;
    [~, class_ind]= max(pred_freq);
    pred_lbl(si,1) = class_list(class_ind);
end
end
   