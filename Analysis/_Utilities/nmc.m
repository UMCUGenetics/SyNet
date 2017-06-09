function [pred, n_optFeat] = nmc(x_tr, l_tr, x_te, opt_info)
if ~exist('opt_info', 'var'), opt_info=struct(); end
Class_list = unique(l_tr);
n_class = numel(Class_list);
n_feat = size(x_tr, 2);

%% Calculate means
mean_tr = zeros(n_class, n_feat);
for ci=1:n_class
    mean_tr(ci,:) = mean(x_tr(l_tr==Class_list(ci), :), 1);
end

%% Prediction
n_optFeat = n_feat;
if isfield(opt_info, 'iCvPar')
	[~, ~, iCvPar] = unique(opt_info.iCvPar, 'Stable');
	n_fold = max(iCvPar);
	in_auc = zeros(n_fold, n_feat);
	top_auc = 0;
	for si=1:n_feat
		for fi=1:n_fold
			iTr = iCvPar~=fi;
			iTe = iCvPar==fi;
			pred = nmc(x_tr(iTr,1:si), l_tr(iTr), x_tr(iTe,1:si));
			in_auc(fi,si) = getAUC(l_tr(iTe), pred, 50);
		end
		if mean(in_auc(:,si))<top_auc
			n_optFeat = si-1;
			fprintf('[%d/%d] Current AUC is: %0.2f%%, Breaking ...\n', si, n_feat, mean(in_auc(:,si))*100);
			break;
		else
			top_auc = mean(in_auc(:,si));
			fprintf('[%d/%d] Current AUC is: %0.2f%%\n', si, n_feat, top_auc*100);
		end
	end
end
dist = pdist2(mean_tr(:,1:n_optFeat), x_te(:,1:n_optFeat));

%% Evaluation
[~, min_ind] = min(dist, [], 1);
pred = Class_list(min_ind);
end