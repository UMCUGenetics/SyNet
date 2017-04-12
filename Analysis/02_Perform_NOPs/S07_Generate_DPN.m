function S07_Generate_DPN()
clc;
clear;

%% Initialization
data_lst = {'ACES' 'META' 'HAIBE'};
sav_path = '../106_Processing_TopPairs/Network_Files/';

%% Loading Top Pairs
for di=1:numel(data_lst)
	
	%% Purify the Pair info
	load(['../106_Processing_TopPairs/Score_Files/PairInfo_' data_lst{di} '.mat'], 'Pair_Info');
	Pair_Info(Pair_Info(:,4)==1 & Pair_Info(:,6)<30, :) = [];
	n_pairs = size(Pair_Info,1);
	[~, sind] = sort(Pair_Info(:, 3), 'descend');
	Pair_Info(sind,12) = 1:n_pairs; % Ranked Assigned Score
	[~, sind] = sort(Pair_Info(:, 11), 'descend');
	Pair_Info(sind,13) = 1:n_pairs; % Ranked Predicted Score
	Pair_Info(:,14) = mean(Pair_Info(:,12:13), 2);
	
	[~, sind] = sort(Pair_Info(:, 14), 'ascend'); % Sort based on combined rank
	Pair_Info = Pair_Info(sind, :);
	
	%% Load gene name
	GeneExpression_Path = getGE_Path(data_lst{di});
	load(GeneExpression_Path, 'Gene_Name');
	n_gene = numel(Gene_Name);
	
	%% Generate Adj Matrix
	ind1 = sub2ind([n_gene n_gene], Pair_Info(:,1), Pair_Info(:,2));
	ind2 = sub2ind([n_gene n_gene], Pair_Info(:,2), Pair_Info(:,1));
	
	Net_Adj = zeros(n_gene, 'single');
	Net_Adj(ind1) = max(Pair_Info(:,14)) - Pair_Info(:,14) + 1;
	Net_Adj(ind2) = max(Pair_Info(:,14)) - Pair_Info(:,14) + 1;
	if any(size(Net_Adj)>n_gene), error(); end
	
	sav_name = [sav_path 'CPN_' data_lst{di} '.mat'];
	fprintf('Saving the results in [%s]\n', sav_name);
	save(sav_name, 'Pair_Info', 'Net_Adj', 'Gene_Name');
end
end
function GeneExpression_Path = getGE_Path(ge_name)
switch ge_name
	case 'ACES'
		GeneExpression_Path = '../../../Dataset/ACES/GE.mat';
	case 'GBM'
		GeneExpression_Path = '../../../Dataset/Cancer_Genomics_Browser/CGB_GBM.mat';
	case 'GBMFl'
		GeneExpression_Path = '../95_Filter_GBM/CGB_GBM_FreqFiltered.mat';
	case 'META'
		GeneExpression_Path = '../94_Pairwise_Evaluation_of_Genes_METABRIC/METABRIC_Filtered.mat';
	case 'HAIBE'
		GeneExpression_Path = '../96_Pairwise_Evaluation_of_Genes_Haibe/Haibe_Filtered.mat';
	otherwise
		error();
end
end