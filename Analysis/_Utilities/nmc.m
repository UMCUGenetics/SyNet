function [pred, n_optFeat] = nmc(x_tr, l_tr, x_te, opt_info)
if ~exist('opt_info', 'var'), opt_info=struct(); end

Class_Labels = unique(l_tr);
n_class = numel(Class_Labels);
n_feat = size(x_tr, 2);
if ~isfield(opt_info, 'MAX_N_Feat'), opt_info.MAX_N_Feat=n_feat; end

%% Calculate means
mean_tr = zeros(n_class, n_feat);
for ci=1:n_class
    mean_tr(ci,:) = mean(x_tr(l_tr==Class_Labels(ci), :), 1);
end

%% Prediction
MAX_N_Feat = min([n_feat opt_info.MAX_N_Feat]);
n_optFeat = n_feat;
if isfield(opt_info, 'iCvPar')
	[~, ~, iCvPar] = unique(opt_info.iCvPar, 'Stable');
	n_fold = max(iCvPar);
	in_auc = zeros(n_fold, MAX_N_Feat);
	for fi=1:n_fold
		iTr = iCvPar~=fi;
		iTe = iCvPar==fi;
		for si=1:MAX_N_Feat
			pred = nmc(x_tr(iTr,1:si), l_tr(iTr), x_tr(iTe,1:si));
			in_auc(fi,si) = getAUC(l_tr(iTe), pred, 50);
		end
		[best_val, best_ind] = max(mean(in_auc(1:fi,:),1));
		fprintf('[%02d/%02d] Fold complete; Best AUC is: %0.2f%% using [%d] features ...\n', fi, n_fold, best_val*100, best_ind);
	end
	[~, n_optFeat] = max(mean(in_auc,1));
end
dist = pdist2(mean_tr(:,1:n_optFeat), x_te(:,1:n_optFeat));

%% Evaluation
[~, min_ind] = min(dist, [], 1);
pred = Class_Labels(min_ind);
end