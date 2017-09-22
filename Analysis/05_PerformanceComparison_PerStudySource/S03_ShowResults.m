clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
res_name = 'Collected_Results.mat';
fprintf('Loading results from [%s] ...\n', res_name);
load(res_name);
n_met = numel(Method_lst);
n_src = numel(Source_lst);
n_rat = numel(Ratio_lst);

%% Plotting
for mi=1:n_met
	close all
	figure('Position', [100 100 1000 700]);
	for si=1:n_src
		for sj=1:n_src
			if si==sj, continue; end
			subplot(n_src, n_src, (si-1)*n_src+sj);
			xTickLabel = {};
			step = 0;
			YLim = [inf -inf];
			for ri=1:n_rat
				step = step + 1;
				val = squeeze(auc_mat(mi,si,sj,ri,:));
				boxplotEx(val, step, {}, struct());
				xTickLabel{step} = sprintf('R=%0.2f', Ratio_lst(ri));
				YLim = [min([YLim(1); val]) max([YLim(2); val])];
			end
			title([Source_lst{si} ' vs. ' Source_lst{sj}], 'FontWeight', 'Bold', 'FontSize', 12);
			set(gca, 'XTick', 1:step, 'XTickLabel', xTickLabel, 'XTickLabelRotation', 40, 'XLim', [0 step+1], ...
				'YLim', [0.50 0.65]);
		end
	end
	annotation('textbox', 'String', Method_lst{mi}, 'Position', [0 1 1 -0.3], ...
		'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Center', ...
		'FontWeight', 'Bold', 'FontSize', 20, 'EdgeColor', 'none');
	ginput
end

