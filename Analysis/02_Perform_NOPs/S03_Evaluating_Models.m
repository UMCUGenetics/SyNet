function [result_name, te_auc] = S03_Evaluating_Models(method_lst, ds_id)
%% ####
% if ismac
% 	fprintf('*** Warning!: Running on debug mode.\n');
% 	method_lst = {'Feral'};
% 	ds_id = '1704161431307296';
% end

%% Initialization
% clc;
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath(genpath('../../../../Useful_Sample_Codes/Mutual_Information_Toolbox'));
dataset_path = './Dataset_Files/';
result_path = './Results_Files/';
opt_info.lasso_opt = {'lassoType', 't', 'CV', 5, 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 0};
opt_info.MAX_N_SUBNET = 500;
if isempty(method_lst), method_lst = {'iPark', 'iChuang', 'iTaylor', 'Feral', 'Single'}; end

%% Get CV info
dataset_list = dir([dataset_path '*_DID-' ds_id '*.mat']);
if numel(dataset_list)~=1, error('Missing or duplicate dataset found.. \n[%s]\n', strjoin({dataset_list.name}, ', ')); end
dataset_name = [dataset_path dataset_list.name];
fprintf('Loading dataset [%s] ...\n', dataset_name);
dataset_info = load(dataset_name);

%% Get method name(s)
n_method = numel(method_lst);
for mi=1:n_method
	method_name = method_lst{mi};
	result_name = sprintf([result_path dataset_list.name(1:end-4) '_MSN-%03d_MTN-' method_name '.mat'], opt_info.MAX_N_SUBNET);
	if exist(result_name, 'file'), error('Results are already computed --> [%s] \n', result_name); end
	
	%% Evaluate model
	fprintf('[%d/%d] Evaluating [%s] model ...\n', mi, n_method, method_name);
	switch method_name
		case 'iPark'
			result = perf_iPark(dataset_info, opt_info);
		case 'iChuang'
			result = perf_iChuang(dataset_info, opt_info);
		case 'iTaylor'
			result = perf_iTaylor(dataset_info, opt_info);
		case 'Feral'
			result = perf_Feral(dataset_info, opt_info);
		case 'Single'
			result = perf_Single(dataset_info, opt_info);
		otherwise
			error('Unknown Method.');
	end
	te_auc = result.te_auc;
	result.Method_Name = method_name;
	result.Dataset_Name = dataset_name;
	
	%% Saving results
	fprintf('Saving results in [%s].\n', result_name);
	save(result_name, '-struct', 'result');
end
end
