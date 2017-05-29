clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop');
addpath('../../../../Useful_Sample_Codes/FisherExactTest/');

%% Load NOPs results
load('../09_Collect_NOPs_Results/Collected_NOP_CV01_Markers.mat', 'Gene_Name', 'Result_Marker', 'method_lst', 'net_lst', 'Result_AUC');
% 1:Method_Index 2:Network_Index 3:Study_Index 4:Repeat_Index 5: Te_AUC
n_met = numel(method_lst);
n_net = numel(net_lst);
n_study = max(Result_Marker(:,3));
n_rep = max(Result_Marker(:,4));
n_gene = numel(Gene_Name);

%% Collecting top genes
Top_List = cell(n_met, n_net, n_study);
for mi=1:n_met
	for ni=1:n_net
		% fprintf('Calculating [%s, %s]\n', method_lst{mi}, net_lst{ni});
		for si=1:n_study
			is_in = Result_Marker(:,1)==mi & Result_Marker(:,2)==ni & Result_Marker(:,3)==si;
			Mrk_set = Result_Marker(is_in,6:end);
			Top_Gene = nonzeros(Mrk_set);
			[Top_List{mi,ni,si}, Top_Freq] = getTop(Top_Gene, 100);
			if numel(Top_List{mi,ni,si})<100
				Top_List{mi,ni,si} = [Top_List{mi,ni,si}; randi(n_gene, 100-numel(Top_List{mi,ni,si}), 1)];
			end
		end
	end
end

%% Measure stability
Study_fet = cell(n_met, n_net);
All_Gene = (1:n_gene)';
for mi=1:n_met
	for ni=1:n_net
		for si=1:n_study
			for sj=si+1:n_study
				[~, ~, pval] = getFET(Top_List{mi,ni,si}, Top_List{mi,ni,sj}, All_Gene);
				Study_fet{mi,ni} = [Study_fet{mi,ni} -log10(pval)];
			end
		end
	end
end

%% Combine P-values
step = 1;
X_lbl = {};
hold on
clr_map = jet(5)*0.9;%[AdvancedColormap('cb', n_net/2); AdvancedColormap('mr', n_net/2)];
for ni=[1 2 4 5 3]
	for mi=2:n_met
		box_h = boxplot(Study_fet{mi,ni}, 'Position', step, 'Color', clr_map(ni,:), 'Widths', 0.5);
		set(box_h, 'LineWidth', 2);
		X_lbl{step,1} = sprintf('%s', method_lst{mi});
		step = step + 1;
	end
	text(step-1, 220, net_lst{ni}, 'HorizontalAlignment', 'Right', 'VerticalAlignment', 'Top', ...
		'FontWeight', 'Bold', 'FontSize', 14);
	step = step + 1;
end
set(gca, 'Xlim', [0 step], 'Ylim', [0 230], ...
	'XTick', 1:step-2, 'XTicklabel', X_lbl, 'XTicklabelRotation', 30, 'FontWeight', 'Bold');
close all

Stab_fet = zeros(n_met, n_net);
for mi=1:n_met
	for ni=1:n_net
		Stab_fet(mi,ni) = median(Study_fet{mi,ni});
	end
end
boxplot(Stab_fet);
set(gca, 'XTick', 1:n_net, 'XTickLabel', net_lst, 'XTickLabelRotation', 25, 'FontWeight', 'Bold', 'FontSize', 14);
ylabel('-log10(pval)', 'FontWeight', 'Bold');

plot(Stab_fet', 'LineWidth', 2);
set(gca, 'XTick', 1:n_net, 'XTickLabel', net_lst, 'XTickLabelRotation', 25, 'FontWeight', 'Bold', 'FontSize', 14);
legend(method_lst);

%% Predictibility
study_auc = median(Result_AUC,4, 'omitnan');
net_auc = median(study_auc, 3, 'omitnan');
boxplot(net_auc);
set(gca, 'XTick', 1:n_net, 'XTickLabel', net_lst, 'XTickLabelRotation', 25, 'FontWeight', 'Bold', 'FontSize', 14);
ylabel('AUC', 'FontWeight', 'Bold');

