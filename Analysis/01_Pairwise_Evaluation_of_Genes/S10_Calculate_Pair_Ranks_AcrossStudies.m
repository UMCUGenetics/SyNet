clc;
clear;

%% Initialization
n_study = 15;
addpath('../../../../Useful_Sample_Codes/getTop');
study_name = cellstr(num2str((1:n_study)','./Network_Files/DSN_SyNetS%02d.mat'));
study_name{n_study} = './Network_Files/DSN_SyNet.mat';
n_top = 10000;

%% Load top pairs
load('./Top_Pairs/Top_SyNet.mat', 'PP_PerStudy', 'PP_Info');
n_pair = size(PP_PerStudy, 1);

%% Main loop
PR_PerStudy = [PP_PerStudy(:,1:2) zeros(n_pair, n_study)];
for si=1:n_study
	fprintf('Loading [%s]...\n', study_name{si});
	load(study_name{si}, 'Net_Adj');
	if ~exist('PP_Index', 'var')
		n_gene = size(Net_Adj,1);
		PP_Index = sub2ind([n_gene n_gene], PP_PerStudy(:,1), PP_PerStudy(:,2));
	end
	fprintf('Sorting ...\n');
	[scr_val, scr_ind] = sort(Net_Adj(:), 'Descend');
	Net_Rnk = zeros(n_gene);
	Net_Rnk(scr_ind) = 1:numel(scr_ind);
	Net_Rnk = Net_Rnk / 2;
	PR_PerStudy(:, si+2) = Net_Rnk(PP_Index);
end

%% Append to 
save('./Top_Pairs/Top_SyNet.mat', 'PR_PerStudy', '-append');