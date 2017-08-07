clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/getRankAUC/');
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/FisherExactTest/');
gs_path = './Enrichment_Sets/';

%% Load NOPs results
load('../S09_Collect_NOPs_Results/Collected_NOP_Markers.mat', 'Gene_Name', 'Result_Marker', 'method_lst', 'net_lst', 'Result_AUC');
% 1:Method_Index 2:Network_Index 3:Study_Index 4:Repeat_Index 5: Te_AUC
n_met = numel(method_lst);
n_net = numel(net_lst);
n_study = max(Result_Marker(:,3));
n_rep = max(Result_Marker(:,4));
n_gene = numel(Gene_Name);

%% Load gene set
gset_file = dir([gs_path 'RefSet_*.csv']);
n_gset = numel(gset_file);
gset_list = cell(n_gset,1);
gset_name = cell(n_gset,1);
for si=1:n_gset
    gset_list{si,1} = regexp(fileread([gs_path gset_file(si).name]), '\n', 'split')';
    file_info = regexp(gset_file(si).name, '[_\.]', 'split');
    gset_name{si,1} = file_info{2};
end

%% Collecting top genes
Top_List = cell(n_met, n_net, n_study);
for mi=1:n_met
	for ni=1:n_net
		% fprintf('Calculating [%s, %s]\n', method_lst{mi}, net_lst{ni});
		for si=1:n_study
			is_in = Result_Marker(:,1)==mi & Result_Marker(:,2)==ni & Result_Marker(:,3)==si;
			Mrk_set = Result_Marker(is_in,6:end);
			Top_List{mi,ni,si} = nonzeros(Mrk_set);
			Top_Gene = nonzeros(Mrk_set);
			%[Top_List{mi,ni,si}, Top_Freq] = getTop(Top_Gene, 100);
% 			if numel(Top_List{mi,ni,si})<100
% 				Top_List{mi,ni,si} = [Top_List{mi,ni,si}; randi(n_gene, 100-numel(Top_List{mi,ni,si}), 1)];
% 			end
		end
	end
end

%% Measure enrichment score
auc_mat = zeros(n_gset, n_met);
step = 1;
for mi=1:n_met
	for ni=[1 2 4 3]
		Mrk_set = unique(vertcat(Top_List{mi,ni,:}), 'Stable');
		%Mrk_set = getTop(vertcat(Top_List{mi,ni,:}), 100);
		rest_set = setdiff(1:n_gene, Mrk_set);
		n_rest = numel(rest_set);
		GSet = [Gene_Name(Mrk_set); Gene_Name(rest_set(randperm(n_rest)))];
		for si=1:n_gset
			auc_mat(si, step) = getRankAUC(gset_list{si}, GSet);
		end
		X_lbl{step,1} = sprintf('%s,%s', method_lst{mi}, net_lst{ni});
		step = step + 1;
	end
end

%% Plotting
close all
figure('Position', [100 100 1600 500]);
bar(auc_mat, 'grouped');
yTick = 0:0.05:1;
yTickLabel = arrayfun(@(x) num2str(x), yTick, 'UniformOutput', false);
set(gca, 'XTick', 1:n_gset, 'XTicklabel', gset_name, 'XTicklabelRotation', 20, 'yTick', yTick, 'yticklabel', yTickLabel, ...
	'ygrid', 'on', 'fontweight', 'bold', 'fontsize', 16);
legend(X_lbl, 'Location', 'BestOutside', 'Orientation', 'Vertical');
% xlim([0 n_gset]);
ylim([0.5 0.65]);
colormap(jet(4));

%% Save plot
% print('-dpdf', '-r300', ['./Plots/Gene_Enrichment.pdf']);

