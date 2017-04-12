clc;

%% Initialization
result_path = './Results_Files/';

%% Collect the results
path_temp = sprintf([result_path 'CID-*.mat']);
res_lst = dir(path_temp);
n_res = numel(res_lst);
fprintf('[%d] result files were found.\n', n_res);

for ri=1:n_res
	res_name = [result_path res_lst(ri).name];
	fprintf('[%03d/%03d] Reading [%s]\n', ri, n_res, res_name);
	
	cv_info = regexp(res_name, '.*CID-(\d+)_CVT(\d+)_ETr-(\w+)_ETe-(\w+)_NET-([\w-.]+)_DID-(\d+)_MSN-100_MTN-(\w+).mat', 'tokens');
	res_data = load(res_name, 'te_auc');
	res_info(ri).CV_ID = cv_info{1}{1};
	res_info(ri).CV_Type = cv_info{1}{2};
	res_info(ri).Tr_Name = cv_info{1}{3};
	res_info(ri).Te_Name = cv_info{1}{4};
	res_info(ri).Net_Name = cv_info{1}{5};
	res_info(ri).Data_ID = cv_info{1}{6};
	res_info(ri).Method_Name = cv_info{1}{7};
	res_info(ri).Data_Name = [res_info(ri).Tr_Name '_' res_info(ri).Te_Name];
	res_info(ri).TE_AUC = res_data.te_auc;
end

%% Save
save('./Collected_Results.mat', 'res_info');
fprintf('All result infos are saved.\n');


