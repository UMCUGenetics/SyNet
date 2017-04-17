clc;

clc;
% clear;

%% Initialization
result_path = './Results_Files/';
net_lst = {'Corr', 'AbsCorr', 'Random', 'String', 'Multinet', 'KEGG', 'MSigDB', ...
	'DSN-SyNet-P99.90', 'DSN-SyNet-P99.99', 'DSN-SyNet-P99.999', ...
	'DSN-SyNet-T00500', 'DSN-SyNet-T01000', 'DSN-SyNet-T05000', 'DSN-SyNet-T10000', 'DSN-SyNet-T20000'
	};
n_net = numel(net_lst);
method_lst = {'iPark', 'iChuang', 'iTaylor', 'Feral', 'Single'};
n_met = numel(method_lst);
cv_id = '170412';
n_study = 14;
n_rep = 20;

te_auc = nan(n_met, n_net, n_study, n_rep);
for si=1:n_study
	for ri=1:n_rep
		for ni=1:n_net
			for mi=1:n_met
				net_name = strrep(net_lst{ni}, 'SyNet', sprintf('SyNetS%02d',si));
				result_temp = [result_path 'CID-' cv_id 'Si' num2str(si,'%02d') 'Ri' num2str(ri,'%02d') '*_NET-' net_name '_*MSN-500_MTN-' method_lst{mi} '.mat'];
				ds_info = dir(result_temp);
				if numel(ds_info)~=1
					fprintf('*** Ignoring [%s]\n', result_temp);
					continue;
				end
				fprintf('Reading from: [%s]\n', ds_info.name);
				data = load([result_path ds_info.name]);
				switch method_lst{mi}
					case 'iPark'
						Dataset_info = load(data.Dataset_Name, 'DatasetTr');
						snet_scr = abs(data.B(:,data.fit.IndexMinMSE));
						n_snet = numel(snet_scr);
						for si=1:n_snet
							data.SubNet_Trimmed{si};
						end
					case 'iChuang'
					case 'iTaylor'
					case 'Feral'
					case 'Single'
					otherwise
						error();
				end
				te_auc(mi,ni,si,ri) = round(data.te_auc*100, 3);
			end
		end
	end
end
save('./Collected_result.mat', 'te_auc');

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


