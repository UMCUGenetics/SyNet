function S05_Collect_PWR(ge_name)
% Run: S05_Collect_PWR('SyNet')

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
pwr_path = ['./PWR_Files/' ge_name '/'];
net_path = './Network_Files/';

%% Load data
pwr_ptrn = [pwr_path 'PWR_' ge_name '_*.mat'];
fprintf('Looking for [%s]\n', pwr_ptrn);
pwr_lst = dir(pwr_ptrn);
n_res = numel(pwr_lst);
fprintf('Found [%d] results ...\n', n_res);

%% Main loop
for fi=1:n_res
	pwr_name = [pwr_path pwr_lst(fi).name];
	fprintf('[%03d/%03d] Reading [%s], ', fi, n_res, pwr_name);
	
	pwr_data = load(pwr_name);
	if ~exist('anc_data', 'var')
		anc_data = pwr_data;
		Gene_Name = anc_data.Gene_Name;
		n_gene = numel(Gene_Name);
		Pair_AUC = zeros(n_gene);
	end
	if ~isequal(anc_data.Patient_Label, pwr_data.Patient_Label) || ...
	   ~isequal(anc_data.Gene_Name, pwr_data.Gene_Name) || ...
	   ~isequal(anc_data.cv_name, pwr_data.cv_name) || ...
	   ~isequal(pwr_data.pair_list, pwr_data.auc_pair(:,1:2))
	   error('Data is not consistent.\n');
	end
	
	auc_pair = [pwr_data.auc_pair(:,1:2) mean(pwr_data.auc_pair(:,3:end),2)];
	n_pair = size(auc_pair, 1);
	fprintf('got [%d] pair AUCs.\n', n_pair);
	for pi=1:n_pair
		Pair_AUC(auc_pair(pi,1), auc_pair(pi,2)) = auc_pair(pi,3);
		Pair_AUC(auc_pair(pi,2), auc_pair(pi,1)) = auc_pair(pi,3);
	end
end
if ~ismac && any(Pair_AUC(:)<0.5)
	fprintf('[i] Warning: Some pairs are missing AUC.\n');
	save('./tmp.mat', 'Pair_AUC');
end
ind_auc = Pair_AUC(1:n_gene+1:end)';

%% Normalization
fprintf('Compute the axis.\n');
pair_max = bsxfun(@max, ind_auc, ind_auc');

x_axis = Pair_AUC;
ox = x_axis-min(x_axis(:));
ox = ox./max(ox(:));

y_axis = Pair_AUC./pair_max;
oy = y_axis-min(y_axis(:));
oy = oy./max(oy(:));

%% Distance calculation
Net_Adj = single(-sqrt((ox-1).^2 + (oy-1).^2));

%% Saving the network
sav_name = [net_path 'DSN_' ge_name '-Mat.mat'];
fprintf('Saving the network in [%s]\n', sav_name);
save(sav_name, 'Net_Adj', 'Pair_AUC', 'Gene_Name');

end