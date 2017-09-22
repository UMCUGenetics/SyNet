clc;
clear;

%% Initialization
res_path = './S01_Performance_PerSource/';
Method_lst = {'Lasso' 'LassoEx'};
n_met = numel(Method_lst);
Source_lst = {'ACES' 'METABRIC' 'TCGA'};
n_src = numel(Source_lst);
Ratio_lst = [0.3 0.5 0.7 0.9 1.0];
n_rat = numel(Ratio_lst);
n_epoch = 5;
auc_mat = nan(n_met, n_src, n_src, n_rat, n_epoch);

%% Main loop
for mi=1:n_met
	for si=1:n_src
		for sj=1:n_src
			if si==sj
				continue;
			end
			for ri=1:n_rat
				res_pattern = sprintf([res_path 'Res_%s_%s_%s_SR%0.2f_*.mat'], Method_lst{mi}, Source_lst{si}, Source_lst{sj}, Ratio_lst(ri));
				res_lst = dir(res_pattern);
				for ei=1:n_epoch
					res_name = [res_path res_lst(ei).name];
					fprintf('Loading [%s] ...\n', res_name);
					res_data = load(res_name, 'te_auc');
					auc_mat(mi, si, sj, ri, ei) = res_data.te_auc;
				end
			end
		end
	end
end

%% Saving
sav_name = 'Collected_Results.mat';
fprintf('Saving results in [%s] ...\n', sav_name);
save(sav_name, 'auc_mat', 'Method_lst', 'Source_lst', 'Ratio_lst');

%{
%% Plotting
boxplot(AUC_mat');
set(gca, 'XTick', 1:11, 'XTickLabel', Method_lst, 'XTickLabelRotation', 40);
plot(AUC_mat(9,:))
%}