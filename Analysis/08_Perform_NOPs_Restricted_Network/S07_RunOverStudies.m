function S07_RunOverStudies(MAX_SUBNET_SIZE, MAX_N_SUBNET)

%% Initialization
if ismac
	MAX_SUBNET_SIZE = 5;
	MAX_N_SUBNET = 10000;
end
clc
data_path = './Dataset_Files/';
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../_Utilities/');

Method_lst = {'RIFeral'};
n_method = numel(Method_lst);
n_study = 12;
opt_info.lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 0};
opt_info.MAX_SUBNET_SIZE = MAX_SUBNET_SIZE;
opt_info.MAX_N_SUBNET = MAX_N_SUBNET;
fprintf('[i] Max SubNet size is: %d\n', MAX_SUBNET_SIZE);
fprintf('[i] Max #SubNet is: %d\n', MAX_N_SUBNET);

%% Main loop
auc_mat = zeros(n_study, n_method);
for si=1:n_study
	data_name = dir([data_path sprintf('DID_Si%02d-Ri001_SyNet-T10000_*.mat', si)]);
	data_name = [data_path data_name.name];
	fprintf('[i] Loading [%s]\n', data_name);
	dataset_info = load(data_name);
	for mi=1:n_method
		switch Method_lst{mi}
			case 'Feral'
				result = perf_Feral(dataset_info, opt_info, 'Mean');
			case 'RIFeral'
				result = perf_Feral(dataset_info, opt_info, 'RI');
			case 'Lasso'
				result = perf_Lasso(dataset_info, opt_info);
		end
		auc_mat(si, mi) = result.te_auc;
	end
end

disp(auc_mat);
disp(mean(auc_mat));