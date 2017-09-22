clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
res_name = 'Collected_Results.mat';
fprintf('Loading results from [%s] ...\n', res_name);
load(res_name);
n_met = numel(Method_lst);
n_std = 12;

%% Collection
auc_index = [];
auc_lst = [];
step = 1;
for si=1:n_std
	for ri=1:5
		auc_index(step, 1:2) = [si ri];
		for mi=1:n_met
			auc_lst(step, mi) = auc_mat(mi, si, ri);
		end
		step = step + 1;
	end
end

%% Comparison
scr_mat = zeros(n_met);
for mi=1:n_met
	for mj=1:n_met
		[h, pval, ci, stats] = ttest(auc_lst(:,mi), auc_lst(:,mj));
		if stats.tstat>0
			scr_mat(mi,mj) = -log10(pval);
		else
			scr_mat(mi,mj) = log10(pval);
		end
	end
end

%% Plotting
close all
img_h = imagesc(scr_mat);
set(img_h, 'AlphaData', ~isnan(scr_mat));
clr_h = colorbar();
clr_map = AdvancedColormap('rwb', 15);
colormap(clr_map);
ylabel(clr_h, '-Log10(pval)');
set(gca, 'XTick', 1:n_met, 'XTickLabel', Method_lst, 'XTickLabelRotation', 35, ...
		 'YTick', 1:n_met, 'YTickLabel', Method_lst, 'YDir', 'Normal', ...
		 'FontWeight', 'Bold', 'FontSize', 12);
caxis([-10 10]);

%% Saving
output_name = sprintf('./Plots/S05_PVal_AcrossStudies.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [36 30], 'PaperPosition', [0 0 36 30]);
print('-dpdf', '-r300', output_name);
fprintf('\n');
%}

%{
for mi=1:n_met
	close all
	X_Label = {};
	step = 1;
	for si=1:n_std
		val = squeeze(auc_mat(mi, si, :));
		val(isnan(val)) = [];
		boxplotEx(val, step, {}, struct());
		X_Label{step, 1} = sprintf('%s', Study_Name{si});
		step = step + 1;
	end
	set(gca, 'XLim', [0 step], 'YLim', [0.50 0.81]);
	set(gca, 'XTick', 1:step-1, 'XTickLabel', X_Label, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', 'FontSize', 12);
	ylabel('AUC', 'FontSize', 12, 'FontWeight', 'Bold');
	title(Method_lst{mi}, 'FontSize', 16, 'FontWeight', 'Bold');
	
	%% Saving
	output_name = sprintf('./Plots/S04_MethodAUC_%s.pdf', Method_lst{mi});
	set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [50 12], 'PaperPosition', [0 0 50 12]);
	print('-dpdf', '-r300', output_name);
	fprintf('\n');
end
%}