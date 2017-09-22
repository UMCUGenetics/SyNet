clc;
clear;

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop');
addpath('../../../../Useful_Sample_Codes/FisherExactTest/');

%% Load NOPs results
load('../02_Perform_NOPs/Collected_NOP_results.mat', 'Gene_Name', 'Result_Info', 'method_lst', 'net_lst', 'Result_AUC');
% 1:Method_Index 2:Network_Index 3:Study_Index 4:Repeat_Index
Result_Info(Result_Info(:,3)>5, :) = [];

n_met = numel(method_lst);
n_net = numel(net_lst);
n_study = max(Result_Info(:,3));
n_rep = max(Result_Info(:,4));
n_gene = numel(Gene_Name);

%% Collecting top genes
Top_List = cell(n_met, n_net, n_study);
for mi=1:n_met
	for ni=1:n_net
		% fprintf('Calculating [%s, %s]\n', method_lst{mi}, net_lst{ni});
		for si=1:n_study
			is_in = Result_Info(:,1)==mi & Result_Info(:,2)==ni & Result_Info(:,3)==si;
			Top_Gene = nonzeros(Result_Info(is_in,6:end));
			[Top_List{mi,ni,si}, Top_Freq] = getTop(Top_Gene, 50);
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
legend(method_lst);

%% Predictibility
study_auc = median(Result_AUC,4, 'omitnan');
net_auc = median(study_auc, 3, 'omitnan');
boxplot(net_auc);
set(gca, 'XTick', 1:n_net, 'XTickLabel', net_lst, 'XTickLabelRotation', 25, 'FontWeight', 'Bold', 'FontSize', 14);
ylabel('AUC', 'FontWeight', 'Bold');

