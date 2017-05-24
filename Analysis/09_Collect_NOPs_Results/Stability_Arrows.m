clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
res_name = './Collected_NOP_results.mat';
fprintf('Loading results from [%s] ...\n', res_name);
load(res_name);
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Selection
mi = 7;
ni = [5 8];
auc_mat = Result_AUC(mi, ni, :, :);
Method_lst = method_lst(mi);
Net_lst = net_lst(ni);
Study_Name = strrep(Study_Name, 'ACES;', '');
n_met = numel(Method_lst);
n_net = numel(Net_lst);
n_std = 12;
Study_Name = Study_Name(1:12);
n_rep = 10;
if any(isnan(auc_mat(:))), error(); end
clr_map = jet(n_std)*0.8;

%% Plotting
for mi=1:n_met
	X_Label = {};
	close all
	hold on
	ln_h = [];
	for si=1:n_std
		for ni=1:n_net
			val_net1 = squeeze(auc_mat(mi, 1, si, :));
			val_net2 = squeeze(auc_mat(mi, 2, si, :));
			
			med_n1 = median(val_net1);
			med_n2 = median(val_net2);
			std_n1 = std(val_net1);
			std_n2 = std(val_net2);
			ln_h(si) = plot([med_n1 med_n2], [std_n1 std_n2], 'k', 'Color', clr_map(si,:), 'LineWidth', 1);
			hold on
			plot(med_n1, std_n1, 'kO', 'MarkerSize', 10, 'Color', clr_map(si,:));
			plot(med_n2, std_n2, 'kx', 'MarkerSize', 15, 'LineWidth', 4, 'Color', clr_map(si,:));
		end
	end
	legend(ln_h, Study_Name, 'FontSize', 12, 'FontWeight', 'Bold', 'Location', 'BestOutside');
	
	set(gca, 'XLim', [0.5 0.8]);
	set(gca, 'FontWeight', 'Bold', 'FontSize', 12);
	xlabel('Median', 'FontSize', 12, 'FontWeight', 'Bold');
	ylabel('Std', 'FontSize', 12, 'FontWeight', 'Bold');
	title({Method_lst{mi} strjoin(Net_lst, ', ')}, 'FontSize', 16, 'FontWeight', 'Bold');
	
	%% Saving
	output_name = sprintf('Arrows_%s_%s-%s.pdf', Method_lst{1}, Net_lst{1}, Net_lst{2});
	set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [25 15], 'PaperPosition', [0 0 25 15]);
	print('-dpdf', '-r300', output_name);
	fprintf('\n');
end
