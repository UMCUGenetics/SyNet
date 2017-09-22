clc;
clear;

%% Initialization
res_path = './Results/';
file_lst = dir([res_path '*.mat']);
n_file = numel(file_lst);
n_study = 14;
Method_lst = {
	'Single-5FCV'
	'Single-LOSO'
	'Single-LOSOMat'
	'SyNetGene-01000-LOSO'
	'SyNetGene-05000-LOSO'
	'SyNetGene-10000-LOSO'
	'SyNetGene-20000-LOSO'
	'SyNetPair-01000-LOSO'
	'SyNetPair-05000-LOSO'
	'SyNetPair-10000-LOSO'
	'SyNetPair-20000-LOSO'
	};
n_met = numel(Method_lst);
AUC_mat = zeros(n_met, n_study);

%% Main loop
for fi=1:n_file
	f_name = [res_path file_lst(fi).name];
	res_info = regexp(f_name, 'Res_Std-(\d+)_Met-([^\.]+)\.mat', 'tokens');
	Study_Ind = str2double(res_info{1}{1});
	Method_Name = res_info{1}{2};
	
	mi = find(ismember(Method_lst, Method_Name));
	if ~isempty(mi)
		if numel(mi)~=1, error(); end
		fprintf('Reading [%s] ... \n', f_name);
		res_data = load(f_name);
		AUC_mat(mi, Study_Ind) = res_data.te_auc;
	end
end

%% Saving
sav_name = 'Collected_Results.mat';
fprintf('Saving results in [%s] ...\n', sav_name);
save(sav_name, 'AUC_mat', 'Method_lst');

%{
%% Plotting
boxplot(AUC_mat');
set(gca, 'XTick', 1:11, 'XTickLabel', Method_lst, 'XTickLabelRotation', 40);
plot(AUC_mat(9,:))
%}