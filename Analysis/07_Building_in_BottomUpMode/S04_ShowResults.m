clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
res_name = 'Collected_Results.mat';
fprintf('Loading results from [%s] ...\n', res_name);
load(res_name);
n_met = numel(Method_lst);
n_std = 12;

%% Plotting
close all
figure('Position', [100 100 1700 400]);
hold on
% %{

for mi=1:n_met
	val = squeeze(auc_mat(mi, :, :));
	boxplotEx(median(val,2, 'omitnan'), mi, {}, struct());
end

set(gca, 'XLim', [0 n_met+1], 'YLim', [0.50 0.75]);
set(gca, 'XTick', 1:mi, 'XTickLabel', Method_lst, 'XTickLabelRotation', 20, 'FontWeight', 'Bold', 'FontSize', 12);
ylabel('AUC', 'FontSize', 12, 'FontWeight', 'Bold');
	
%% Saving
output_name = sprintf('./Plots/S04_OverallAUC_AcrossStudies.pdf');
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [60 15], 'PaperPosition', [0 0 60 15]);
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