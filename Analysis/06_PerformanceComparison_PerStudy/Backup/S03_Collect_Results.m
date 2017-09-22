clc;
clear
close all

Study_lst = {'ACES' 'METABRIC' 'TCGA'};
auc_mat = nan(3);
for si=1:3
	for sj=1:3
		if si==sj, continue; end
		dt_name = sprintf('./S01_Matlab_Study_Test/Lasso_%s_%s.mat', Study_lst{si}, Study_lst{sj});
		res_data = load(dt_name, 'te_auc');
		auc_mat(si,sj) = res_data.te_auc;
	end
end

imagesc(auc_mat);
% caxis([0.5 0.7]);
clr_h = colorbar();
xlabel('Testing');
ylabel('Training');
ylabel(clr_h, 'AUC');
set(gca, 'XTick', 1:3, 'XTickLabel', Study_lst, 'XTickLabelRotation', 45, ...
	'YTick', 1:3, 'YTickLabel', Study_lst, 'FontWeight', 'Bold');