clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
Study_Name = strrep(Study_Name, 'ACES;', '');

%% Loading results
Method_lst = {};
Net_lst = {};
AUC_cell = {};
step = 1;
for cv_id=[91 1]
	res_name = sprintf('Collected_NOP_AUCs_CV%02d.mat', cv_id);
	fprintf('Loading results from [%s] ...\n', res_name);
	clc_data = load(res_name);
	for ni=1:numel(clc_data.net_lst)
		for si=1:size(clc_data.Result_AUC,3)
			for mi=1:numel(clc_data.method_lst)
				AUC_cell{step,1} = squeeze(clc_data.Result_AUC(mi,ni,si,:));
				Met_lst{step,1} = sprintf('%d-%s-%s-%d', cv_id, clc_data.method_lst{mi}, clc_data.net_lst{ni}, si);
				Met_Src{step,1} = Study_Name{si};
				step = step + 1;
			end
		end
	end
end

%% Select Methods
net_name = clc_data.net_lst{2};
% ptr_str = sprintf('.*-%s-%s-.*', 'T(Ag|)NMC', net_name); 
% ptr_str = sprintf('.*-%s-%s-.*', 'TAg(NMC|LEx)', net_name);
ptr_str = sprintf('.*-%s-%s-.*', 'TAg(NMC)', net_name);
% ptr_str = sprintf('91-TAgNMC.*');
is_in = cellfun('length', regexp(Met_lst, ptr_str));
Met_lst = Met_lst(is_in>0);
AUC_cell = AUC_cell(is_in>0);
Met_Src = Met_Src(is_in>0);
n_met = numel(Met_lst);

%% Plot
close all
figure('Position', [100 100 1500 400]);
hold on
for mi=1:n_met
	switch Met_lst{mi}(end) %mod(mi-1,4)+1
	%switch num2str(mi)
		case cellstr(num2str([1 5]'))
			Point_Param.Colormap = AdvancedColormap('cb', 10)*0.8;
		case cellstr(num2str([2 6 9]'))
			Point_Param.Colormap = AdvancedColormap('lg', 10)*0.8;
		case cellstr(num2str([3 7]'))
			Point_Param.Colormap = AdvancedColormap('mr', 10)*0.8;
		case cellstr(num2str([4 8 0]'))
			Point_Param.Colormap = AdvancedColormap('yo', 10)*0.8;
	end
	Point_Param.ColorCAxis = [min(AUC_cell{mi}) max(AUC_cell{mi})];
	box_h = boxplotEx(AUC_cell{mi}, mi, {}, Point_Param);
	set(box_h, 'Color', Point_Param.Colormap(5,:), 'Marker', 'none');
end
set(gca, 'XLim', [0 n_met+1], 'YLim', [0.50 0.75]);
set(gca, 'XTick', 1:n_met, 'XTickLabel', Met_Src, 'XTickLabelRotation', 15, 'FontWeight', 'Bold', 'FontSize', 10);
ylabel('AUC', 'FontSize', 14, 'FontWeight', 'Bold');
title(ptr_str, 'FontSize', 14);

%% Saving
output_name = sprintf('./Plots/S04_AUC_Stability_Net-%s_%s.pdf', net_name, ptr_str);
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [35 12], 'PaperPosition', [0 0 35 12]);
print('-dpdf', '-r300', output_name);
fprintf('\n');
