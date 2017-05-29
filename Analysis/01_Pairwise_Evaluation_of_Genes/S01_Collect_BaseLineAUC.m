clc;
clear;

%% Initialization
result_path = '../08_Perform_NOPs_Restricted_Network/Results_Files/';
Method_Name = 'TAgNMC';
Net_Name = 'Random-T00020';
cv_ind = 1;
n_study = 14;
n_rep = 10;
sav_path = './Baseline_AUCs/';
[~,~] = mkdir(sav_path);

%% Main loop
NMC_AUC = nan(n_study, n_rep);
for si=1:n_study
	for ri=1:n_rep
		res_name = sprintf('%sDID_CVT%02d_Si%02d-Ri%03d_%s_*_MTN-%s.mat', result_path, cv_ind, si, ri, Net_Name, Method_Name);
		file_info = dir(res_name);
		if numel(file_info)~=1, error(); end
		fprintf('Reading from: [%s]\n', file_info.name);
		res_data = load([result_path file_info.name]);
		NMC_AUC(si, ri) = res_data.te_auc;
	end
end
if any(isnan(NMC_AUC(:)))
	fprintf('[i] Warning some values are missing.\n');
end
BaseLine_AUC = median(NMC_AUC, 2);

%% Load CV ind
load('../01_Pairwise_Evaluation_of_Genes/CV_Files/CV_SyNet_CVT02.mat', 'cv_obj');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Index');
n_rep = size(cv_obj, 2);
CVB_AUC = zeros(n_study, n_rep);
for si=1:n_study
	for ri=1:n_rep
		Te_ind = unique(Study_Index(cv_obj(si,ri).iTe));
		if numel(Te_ind)~=1, error(); end
		CVB_AUC(si, ri) = BaseLine_AUC(Te_ind);
	end
end

%% Saving
out_name = sprintf([sav_path 'BA_CV%02d_%s_%s.mat'], cv_ind, Method_Name, Net_Name);
save(out_name, 'BaseLine_AUC', 'CVB_AUC', 'Method_Name', 'Net_Name');
