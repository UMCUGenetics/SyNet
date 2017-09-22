clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
res_name = './Collected_NOP_CV01_AUCs.mat';
fprintf('Loading results from [%s] ...\n', res_name);
load(res_name);
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
n_std = 12;
Study_Name = strrep(Study_Name(1:n_std), 'ACES;', '');
n_rep = 10;
mi = 5;

%% Load base line
load('../01_Pairwise_Evaluation_of_Genes/Baseline_AUCs/BA_CV01_TAgNMC_Random-T00020.mat');

%% Selection
clr_map = jet(n_std)*0.8;
Method_Name = method_lst{mi};

%% Plotting
for ni=[12 14 15 16]
	net_ind = [ni 7];
	auc_mat = Result_AUC(mi, net_ind, :, :);
	if any(isnan(auc_mat(:))), error(); end
	X_Label = {};
	close all
	hold on
	ln_h = [];
	for si=1:n_std
		val_net1 = squeeze(auc_mat(1, 1, si, :));
		val_net2 = squeeze(auc_mat(1, 2, si, :));
		
		med_n1 = median(val_net1);
		med_n2 = median(val_net2);
		std_n1 = std(val_net1);
		std_n2 = std(val_net2);
		ln_h(si) = plot([med_n1 med_n2], [std_n1 std_n2], 'k', 'Color', clr_map(si,:), 'LineWidth', 2);
		hold on
		plot(med_n1, std_n1, 'kO', 'MarkerSize', 10, 'Color', clr_map(si,:));
		plot(med_n2, std_n2, 'kx', 'MarkerSize', 15, 'LineWidth', 4, 'Color', clr_map(si,:));
	end
	legend(ln_h, Study_Name, 'FontSize', 12, 'FontWeight', 'Bold', 'Location', 'BestOutside');
	
	set(gca, 'XLim', [0.5 0.8]);
	set(gca, 'FontWeight', 'Bold', 'FontSize', 12);
	xlabel('Median', 'FontSize', 12, 'FontWeight', 'Bold');
	ylabel('Std', 'FontSize', 12, 'FontWeight', 'Bold');
	title({Method_Name strjoin(net_lst(net_ind), ', ')}, 'FontSize', 16, 'FontWeight', 'Bold');
	
	%% Saving
	output_name = sprintf('./Plots/Arrows_%s_%s-%s.pdf', Method_Name, net_lst{net_ind(1)}, net_lst{net_ind(2)});
	set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [25 15], 'PaperPosition', [0 0 25 15]);
	print('-dpdf', '-r300', output_name);
	fprintf('\n');
end

