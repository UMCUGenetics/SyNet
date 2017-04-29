clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
res_name = 'Collected_Results.mat';
fprintf('Loading results from [%s] ...\n', res_name);
load(res_name);
n_src = numel(Source_lst);
n_std = numel(Study_Name);

%% Plotting
close all
figure('Position', [100 100 1500 500]);
for si=1:n_src
	subplot(n_src, 1, si);
	hold on
	for sj=1:n_std
		val = auc_cell{si, sj};
		if isempty(val), continue; end
		boxplotEx(val, sj, {}, struct());
	end
	set(gca, 'XTick', [], 'XLim', [0 n_std+1], 'YLim', [0.50 0.81]);
	text(0, 0.81, Source_lst{si}, 'FontWeight', 'Bold', 'FontSize', 14, 'VerticalAlignment', 'Top', 'HorizontalAlignment', 'Left');
end
set(gca, 'XTick', 1:n_std, 'XTickLabel', Study_Name, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', 'FontSize', 12);


