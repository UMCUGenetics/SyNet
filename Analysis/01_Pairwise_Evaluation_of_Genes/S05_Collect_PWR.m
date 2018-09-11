function S05_Collect_PWR(ge_name)
% Run: S05_Collect_PWR('SyNet')
% sinter bigmem --mem=200GB /opt/insy/matlab/R2018a/bin/matlab

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
if ismac, ge_name='SyNet'; end
pwr_path = ['./PWR_Files/' ge_name '/'];
net_path = './Network_Files/';

%% Load data
pwr_ptrn = [pwr_path 'PWR_' ge_name '_*.mat'];
fprintf('Looking for [%s]\n', pwr_ptrn);
pwr_lst = dir(pwr_ptrn);
n_file = numel(pwr_lst);
fprintf('Found [%d] results ...\n', n_file);

%% Main loop
for fi=1:n_file
	pwr_name = [pwr_path pwr_lst(fi).name];
	fprintf('[%03d/%03d] Reading [%s], ', fi, n_file, pwr_name);
	
	pwr_data = load(pwr_name);
	if ~exist('anc_data', 'var')
		anc_data = pwr_data;
		Gene_Name = anc_data.Gene_Name;
		n_gene = numel(Gene_Name);
		n_study = size(anc_data.auc_pair,2) - 2;
		fprintf('\n[%d] studies are found ...\n', n_study);
		for si=1:n_study+1
			Total_AUC{si,1} = zeros(n_gene);
			Total_Std{si,1} = zeros(n_gene);
		end
		%bl_data = load('./Baseline_AUCs/BA_CV01_TAgNMC_Random-T00020.mat');
	end
	if ~isequal(anc_data.Patient_Label, pwr_data.Patient_Label) || ...
			~isequal(anc_data.Gene_Name, pwr_data.Gene_Name) || ...
			~isequal(anc_data.cv_name, pwr_data.cv_name) || ...
			~isequal(pwr_data.pair_list, pwr_data.auc_pair(:,1:2))
		error('Data is not consistent.\n');
	end
	
	%% Get Std
	auc_cell = pwr_data.auc_cell;
	n_pair = size(auc_cell, 1);
	fprintf('got [%d] pair AUCs.\n', n_pair);
	std_pair = zeros(n_pair, n_study+1);
	for pi=1:n_pair
		std_pair(pi,1:n_study) = std(double(auc_cell{pi})/10000, 0, 2);
	end
	
	%% Get Median
	std_pair(:, end) =               std(pwr_data.auc_pair(:,3:end),0,2);
	auc_pair = [pwr_data.auc_pair median(pwr_data.auc_pair(:,3:end), 2)];
	for si=1:n_study+1
		for pi=1:n_pair
			Total_AUC{si}(auc_pair(pi,1), auc_pair(pi,2)) = auc_pair(pi,si+2);
			Total_AUC{si}(auc_pair(pi,2), auc_pair(pi,1)) = auc_pair(pi,si+2);
			Total_Std{si}(auc_pair(pi,1), auc_pair(pi,2)) = std_pair(pi,si);
			Total_Std{si}(auc_pair(pi,2), auc_pair(pi,1)) = std_pair(pi,si);
		end
	end
end

%% Saving the networks
for si=1:n_study+1
	Pair_AUC = Total_AUC{si};
	Pair_Std = Total_Std{si};
	if ~ismac && any(Pair_AUC(:)<0.5)
		fprintf('[i] Warning: Some pairs are missing AUC in [%d].\n', si);
		save(sprintf('./tmp_%d.mat', si), 'Pair_AUC');
	end
	
	%% Normalization
	fprintf('Compute the axis for [%d] study.\n', si);
	
	x_axis = Pair_AUC;
	ox = x_axis-min(x_axis(:));
	ox = ox./max(ox(:));
	
	ind_auc = Pair_AUC(1:n_gene+1:end)';
	pair_max = bsxfun(@max, ind_auc, ind_auc');
	y_axis = Pair_AUC./pair_max;
	oy = y_axis-min(y_axis(:));
	oy = oy./max(oy(:));
	
	%% Distance calculation
	Net_Adj = single(-sqrt((ox-1).^2 + (oy-1).^2));
	
	%% Saving the network
	sav_name = [net_path 'DSN_' ge_name 'S' num2str(si, '%02d') '.mat'];
	fprintf('Saving network in [%s] ... \n', sav_name);
	save(sav_name, 'Net_Adj', 'Pair_AUC', 'Pair_Std', 'Gene_Name', 'anc_data');
end
fprintf('Revome the numbering from last file. e.g. DSN_SyNet.mat');