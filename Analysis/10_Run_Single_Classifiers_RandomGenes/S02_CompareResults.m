clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
res_name = 'Collected_NOP_AUCs.mat';
fprintf('Loading results from [%s] ...\n', res_name);
clc_data = load(res_name);
Method_lst = clc_data.method_lst;
Net_lst = clc_data.net_lst;
auc_mat = clc_data.Result_AUC;
n_met = numel(Method_lst);
n_net = numel(Net_lst);
n_std = 12;
n_rep = 10;

%{%
auc_mat = mean(auc_mat, 3, 'omitnan');
imagesc(auc_mat);
colormap(jet(10));
% colormap(copper(10));
set(gca, 'XTick', 1:n_net, 'XTickLabel', Net_lst, 'XTickLabelRotation', 40, ...
	'YTick', 1:n_met, 'YTickLabel', Method_lst);
colorbar();
% caxis([0.63 0.68]);
disp;
%}

%{%
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(10);%[AdvancedColormap('cb', n_net/2); AdvancedColormap('mr', n_net/2)];
X_lbl = {};
step = 1;
% auc_mat = median(auc_mat, 4, 'omitnan');
anc_auc = squeeze(auc_mat(2, 5, :, :));
for ni=[4 7 9:10 8]
	for mi=[1 3 5 12]
		auc = squeeze(auc_mat(mi, ni, :, :)); - anc_auc; % 
		[muhat(1), sigmahat, muci, sigmaci] = normfit(repmat(auc(:),10,1));
		muhat(2) = muhat(1) - muci(1);
		muhat(3) = muci(2) - muhat(1);
		
		bar(step, muhat(1), 'FaceColor', clr_map(ni,:));
		errorbarEx(step, muhat(1), muhat(2), muhat(3), 2, 0.2, [0 0 0]);
		X_lbl{step, 1} = sprintf('%s', Method_lst{mi});
		step = step + 1;
	end
	text(step, 0.69, Net_lst{ni}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top', ...
		'FontWeight', 'Bold', 'FontSize', 14);
	step = step + 3;
end
xlim([0 step]);
set(gca, 'XTick', 1:step-4, 'XTickLabel', X_lbl, 'XTickLabelRotation', 20, 'FontWeight', 'Bold');
colormap(clr_map);
ylim([0.61 0.69]);
disp;
% legend(Net_lst, 'FontWeight', 'Bold');
%}

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
		set(gca, 'XLim', [0 step], 'YLim', [50 81]);
		set(gca, 'XTick', 1:step-1, 'XTickLabel', X_Label, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', 'FontSize', 12);
		ylabel('AUC', 'FontSize', 12, 'FontWeight', 'Bold');
		title([Method_lst{mi} ', ' Net_lst{ni}], 'FontSize', 16, 'FontWeight', 'Bold');
		
		%% Saving
		output_name = sprintf('./Plots/S02_AUCDist_Per-Method-Net_%s_%s.pdf', Method_lst{mi}, Net_lst{ni});
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
auc_cell = {};
% auc_mat = mean(auc_mat, 4);
for mi=1:n_met
	for ni=1:n_net
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
clr_map = AdvancedColormap('rwb', 15);
colormap(clr_map);
ylabel(clr_h, '-Log10(pval)');
set(gca, 'XTick', 1:n_elm, 'XTickLabel', X_lbl, 'XTickLabelRotation', 35, ...
		'YTick', 1:n_elm, 'YTickLabel', X_lbl, 'YDir', 'Normal', ...
	'FontWeight', 'Bold', 'FontSize', 7);
caxis([-15 15]);

%% Saving
output_name = sprintf('./Plots/S02_MetNetAUC.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [30 25], 'PaperPosition', [0 0 30 25]);
print('-dpdf', '-r300', output_name);
fprintf('\n');
%}

%{		
for si=1:n_std
	close all
	X_Label = {};
	step = 1;
	for mi=1:n_met
		for ni=1:n_net
			val = squeeze(auc_mat(mi, ni, si, :));
			val(isnan(val)) = [];
			boxplotEx(val, step, {}, struct());
			X_Label{step, 1} = sprintf('%s,%s', Method_lst{mi}, Net_lst{ni});
			step = step + 1;
		end
	end
	set(gca, 'XLim', [0 step], 'YLim', [50 81]);
	set(gca, 'XTick', 1:step-1, 'XTickLabel', X_Label, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', 'FontSize', 12);
	ylabel('AUC', 'FontSize', 12, 'FontWeight', 'Bold');
	title(Study_Name{si}, 'FontSize', 16, 'FontWeight', 'Bold');
	
	%% Saving
	output_name = sprintf('./Plots/S02_AUCDist_Per-Study_%s.pdf', Study_Name{si});
	set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [70 15], 'PaperPosition', [0 0 70 15]);
	print('-dpdf', '-r300', output_name);
	fprintf('\n');
end
%}

%{%
opt_lst = [
	4 6
	2 6
	7 4
	8 5
	];
close all
figure('Position', [100 100 1500 400]);
hold on
X_Label = {};
step = 1;
for si=1:n_std
	val_lst = {};
	for oi=1:size(opt_lst)
		mi = opt_lst(oi,1);
		ni = opt_lst(oi,2);
		
		val = squeeze(auc_mat(mi, ni, si, :));
		val(isnan(val)) = [];
		switch oi
			case 1
				Point_Param.Colormap = AdvancedColormap('cb', 10)*0.8;
			case 2
				Point_Param.Colormap = AdvancedColormap('lg', 10)*0.8;
			case 3
				Point_Param.Colormap = AdvancedColormap('mr', 10)*0.8;
			case 4
				Point_Param.Colormap = AdvancedColormap('yo', 10)*0.8;
		end
		Point_Param.ColorCAxis = [min(val) max(val)];
		box_h = boxplotEx(val, step, {}, Point_Param);
		set(box_h, 'Color', Point_Param.Colormap(5,:), 'Marker', 'none');
		X_Label{step, 1} = sprintf('%s,%s', Method_lst{mi}, Net_lst{ni});
		step = step + 1;
		val_lst{oi} = val;
		if oi>1
			[~, pval] = ttest(val_lst{oi-1}, val_lst{oi});
			text(step-1.5, 80, sprintf('%0.0e', pval), 'FontSize', 14, 'FontWeight', 'Bold', 'Rotation', 90, 'HorizontalAlignment', 'Right');
		end
	end
end
set(gca, 'XLim', [0 step], 'YLim', [0.50 0.81]);
set(gca, 'XTick', 1:step-1, 'XTickLabel', X_Label, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', 'FontSize', 8);
ylabel('AUC', 'FontSize', 12, 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S02_ResDist.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [70 12], 'PaperPosition', [0 0 70 12]);
print('-dpdf', '-r300', output_name);
fprintf('\n');

%}