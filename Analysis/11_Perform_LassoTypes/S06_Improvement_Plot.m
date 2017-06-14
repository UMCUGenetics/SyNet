clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
res_name = 'Collected_Results_MSN-500.mat';
fprintf('Loading results from [%s] ...\n', res_name);
load(res_name);
Study_Name = strrep(Study_Name, 'ACES;', '');
n_met = numel(Method_lst);
n_net = numel(Net_lst);
n_std = 12;
n_rep = 10;

%% Selection
lst = [
	1 2
	1 4
];
auc_cell = {};

%% Generate AUC set
auc_cell = cell(n_std, 1);
for li=1:size(lst, 1)
	Method_Name = Method_lst{lst(li,1)};
	Net_Name = Net_lst{lst(li,2)};
	
	for ni=1:n_net
		val = squeeze(auc_mat(lst(li,1), lst(li,2), si, :));
		val(isnan(val)) = [];
		
		auc_cell{ni, 1} = val;
		box_h = boxplotEx(val, step, {}, Point_Param);
		set(box_h, 'Color', Point_Param.Colormap(5,:), 'Marker', 'none');
		X_Label{step, 1} = sprintf('%s', Study_Name{si});
		step = step + 1;
	end
end
%% Plotting
X_Label = {};
close all
step = 1;
for si=1:n_std
	for ni=1:n_net
		val = squeeze(auc_mat(1, ni, si, :));
		val(isnan(val)) = [];
		switch ni
			case 1
				Point_Param.Colormap = AdvancedColormap('cb', 10)*0.8;
			case 2
				Point_Param.Colormap = AdvancedColormap('lg', 10)*0.8;
			case 3
				Point_Param.Colormap = AdvancedColormap('mr', 10)*0.8;
		end
		Point_Param.ColorCAxis = [min(val) max(val)];
		box_h = boxplotEx(val, step, {}, Point_Param);
		set(box_h, 'Color', Point_Param.Colormap(5,:), 'Marker', 'none');
		X_Label{step, 1} = sprintf('%s', Study_Name{si});
		step = step + 1;
	end
end
set(gca, 'XLim', [0 step], 'YLim', [0.50 0.81]);
set(gca, 'XTick', 1:step-1, 'XTickLabel', X_Label, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', 'FontSize', 12);
ylabel('AUC', 'FontSize', 12, 'FontWeight', 'Bold');
%legend({'Corr' 'SyNet'});
title({Method_lst{1} strjoin(Net_lst, ', ')}, 'FontSize', 16, 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/Stability_%s_%s.pdf', Method_lst{1}, trjoin(Net_lst, '-'));
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [60 17], 'PaperPosition', [0 0 60 17]);
print('-dpdf', '-r300', output_name);
fprintf('\n');

%% Collection
%{
auc_index = [];
auc_lst = [];
step = 1;
for ni=1:n_net
	for si=1:n_std
		for ri=1:n_rep
			auc_index(step, 1:3) = [ni si ri];
			for mi=1:n_met
				auc_lst(step, mi) = auc_mat(mi, ni, si, ri);
			end
			step = step + 1;
		end
	end
end
del_ind = any(isnan(auc_lst),2);
auc_index(del_ind, :) = [];
auc_lst(del_ind, :) = [];
%}

%% Comparison
%{
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
%}

%% Plotting
%{
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
	for ni=1:n_net
		close all
		X_Label = {};
		step = 1;
		for si=1:n_std
			val = squeeze(auc_mat(mi, ni, si, :));
			val(isnan(val)) = [];
			boxplotEx(val, step, {}, struct());
			X_Label{step, 1} = sprintf('%s', Study_Name{si});
			step = step + 1;
		end
		set(gca, 'XLim', [0 step], 'YLim', [0.50 0.81]);
		set(gca, 'XTick', 1:step-1, 'XTickLabel', X_Label, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', 'FontSize', 12);
		ylabel('AUC', 'FontSize', 12, 'FontWeight', 'Bold');
		title([Method_lst{mi} ', ' Net_lst{ni}], 'FontSize', 16, 'FontWeight', 'Bold');
		
		%% Saving
		output_name = sprintf('./Plots/S05_MethodAUC_%s_%s.pdf', Method_lst{mi}, Net_lst{ni});
		set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [50 12], 'PaperPosition', [0 0 50 12]);
		print('-dpdf', '-r300', output_name);
		fprintf('\n');
	end
end
%}

%{
for mi=1:n_met
	scr_mat = zeros(n_net);
	for ni=1:n_net
		for nj=1:n_net
			val_i = squeeze(auc_mat(mi, ni, :, :));
			val_j = squeeze(auc_mat(mi, nj, :, :));
			
			[h, pval, ci, stats] = ttest(val_i(:), val_j(:));
			if stats.tstat>0
				scr_mat(ni,nj) = -log10(pval);
			else
				scr_mat(ni,nj) = log10(pval);
			end
		end
	end
	
	%% Plotting
	close all
	img_h = imagesc(scr_mat);
	set(img_h, 'AlphaData', ~isnan(scr_mat));
	clr_h = colorbar();
	clr_map = AdvancedColormap('rwb', 9);
	colormap(clr_map);
	ylabel(clr_h, '-Log10(pval)');
	title(Method_lst{mi}, 'FontSize', 16, 'FontWeight', 'Bold');
	set(gca, 'XTick', 1:n_net, 'XTickLabel', Net_lst, 'XTickLabelRotation', 35, ...
		'YTick', 1:n_net, 'YTickLabel', Net_lst, 'YDir', 'Normal', ...
		'FontWeight', 'Bold', 'FontSize', 12);
	caxis([-4 4]);
	
	%% Saving
	output_name = sprintf('./Plots/S05_NetAUC_%s_%s.pdf', Method_lst{mi});
	set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [18 14], 'PaperPosition', [0 0 18 14]);
	print('-dpdf', '-r300', output_name);
	fprintf('\n');
end
%}

%{
step = 1;
X_lbl = {};
auc_mat = median(auc_mat, 4);
for ni=1:n_net
	for mi=1:n_met
		auc_cell{step,1} = squeeze(auc_mat(mi, ni, :));
		X_lbl{step, 1} = sprintf('%s,%s', Method_lst{mi}, Net_lst{ni});
		step = step + 1;
	end
end
n_elm = numel(auc_cell);

scr_mat = zeros(n_elm);
for ei=1:n_elm
	for ej=1:n_elm
		
		[h, pval, ci, stats] = ttest(auc_cell{ei}, auc_cell{ej});
		if stats.tstat>0
			scr_mat(ei,ej) = -log10(pval);
		else
			scr_mat(ei,ej) = log10(pval);
		end
	end
end

%% Plotting
close all
img_h = imagesc(scr_mat);
set(img_h, 'AlphaData', ~isnan(scr_mat));
clr_h = colorbar();
clr_map = AdvancedColormap('rwb', 9);
colormap(clr_map);
ylabel(clr_h, '-Log10(pval)');
set(gca, 'XTick', 1:n_elm, 'XTickLabel', X_lbl, 'XTickLabelRotation', 35, ...
		'YTick', 1:n_elm, 'YTickLabel', X_lbl, 'YDir', 'Normal', ...
	'FontWeight', 'Bold', 'FontSize', 12);
caxis([-4 4]);

%% Saving
output_name = sprintf('./Plots/S05_MetNetAUC.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [18 14], 'PaperPosition', [0 0 18 14]);
print('-dpdf', '-r300', output_name);
fprintf('\n');
%}

