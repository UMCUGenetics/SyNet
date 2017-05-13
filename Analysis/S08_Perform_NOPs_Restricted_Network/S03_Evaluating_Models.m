function [result_name, te_auc] = S03_Evaluating_Models(method_lst, ds_id)
%% ####
if ismac && ~exist('method_lst', 'var')
	fprintf('*** Warning!: Running on debug mode.\n');
	method_lst = {'RegAG'};
	ds_id = 'Si01-Ri001_SyNet-T10000';
end

%% Initialization
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../_Utilities/');
dataset_path = './Dataset_Files/';
result_path = './Results_Files/';
opt_info.lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 0};
opt_info.MAX_N_SUBNET = 500;
if isempty(method_lst), method_lst = {'iPark', 'iChuang', 'iTaylor', 'Feral', 'Lasso'}; end

%% Get CV info
dataset_list = dir([dataset_path 'DID_' ds_id '*.mat']);
if numel(dataset_list)~=1, error('Missing or duplicate dataset found.. \n[%s]\n', strjoin({dataset_list.name}, ', ')); end
dataset_name = [dataset_path dataset_list(1).name];
fprintf('Loading dataset [%s] ...\n', dataset_name);
dataset_info = load(dataset_name);

%% Get method name(s)
n_method = numel(method_lst);
te_auc = zeros(n_method,1);
for mi=1:n_method
	method_name = method_lst{mi};
	result_name = sprintf([result_path '%s_MSN-%03d_MTN-%s.mat'], dataset_list(1).name(1:end-4), opt_info.MAX_N_SUBNET, method_name);
	if exist(result_name, 'file')
		fprintf('[i] Results are already computed --> [%s] \n', result_name);
		continue;
	end
	
	%% Evaluate model
	fprintf('[%d/%d] Evaluating [%s] model ...\n', mi, n_method, method_name);
	switch method_name
		case 'iPark'
			result = perf_iPark(dataset_info, opt_info, 'Mean');
		case 'RI-iPark'
			result = perf_iPark(dataset_info, opt_info, 'RI');
		case 'iChuang'
			result = perf_iChuang(dataset_info, opt_info);
		case 'iTaylor'
			result = perf_iTaylor(dataset_info, opt_info);
		case 'Feral'
			result = perf_Feral(dataset_info, opt_info, 'Mean');
		case 'RI-Feral'
			result = perf_Feral(dataset_info, opt_info, 'RI');
		case 'AS-Feral'
			highSN_info = opt_info;
			highSN_info.MAX_N_SUBNET = 10000;
			result = perf_Feral(dataset_info, highSN_info, 'RI');
		case 'Lasso'
			result = perf_Lasso(dataset_info, opt_info);
		case 'LExAG'
			result = perf_LExAG(dataset_info, opt_info);
		case 'TLEx'
			result = perf_TLEx(dataset_info, opt_info);
		case 'TNMC'
			result = perf_TTNMC(dataset_info, opt_info);
		case 'Regress'
			result = perf_Regress(dataset_info, opt_info);
		case 'RegAG'
			result = perf_RegAG(dataset_info, opt_info);
		otherwise
			error('Unknown Method.');
	end
	te_auc(mi) = result.te_auc;
	result.Method_Name = method_name;
	result.Dataset_Name = dataset_name;
	
	%% Saving results
	fprintf('Saving results in [%s].\n', result_name);
	save(result_name, '-struct', 'result');
end
end
