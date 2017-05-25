clc;
clear;

%% Initialization
result_path = '../S08_Perform_NOPs_Restricted_Network/Results_Files/';
method_lst = {'TNMC' 'TAgLEx' 'TLEx' 'LExAG' 'Lasso' 'iPark' 'RI-iPark' 'Feral'  'RI-Feral' 'AS-Feral' 'GLasso' 'CFGLasso'}; %
n_met = numel(method_lst);
net_lst = {
	'SyNet-T10000', 'MinSyn-T10000', 'AvgSyn-T10000', 'CrSyn-T10000', 'CrMinSyn-T10000' ...
	'Corr-T10000','AbsCorr-T10000','Random-T10000', ...
	'STRING-T10000', 'KEGG-T10000'
	};
n_net = numel(net_lst);
n_study = 12;
n_rep = 10;
MAX_N_SN = 500;
sav_path = './Collected_AUCs/';
[~,~] = mkdir(sav_path);

%% Main loop
Result_AUC = nan(n_met, n_net, n_study, n_rep);
for mi=1:n_met
	for ni=1:n_net
		Method_Name = method_lst{mi};
		Net_Name = net_lst{ni};
		sav_name = sprintf([sav_path 'CRES_CVT01_%s_%s_MSN-%03d.mat'], Method_Name, Net_Name, MAX_N_SN);
		if exist(sav_name, 'file')
			fprintf('/// Loading from [%s]\n', sav_name);
			res_data = load(sav_name);
			if ~isequal(Method_Name, res_data.Method_Name) || ~isequal(Net_Name, res_data.Net_Name)|| ~isequal(MAX_N_SN, res_data.MAX_N_SN), error(); end
			Part_AUC = res_data.Part_AUC;
			if any(isnan(Part_AUC(:))), fprintf('Warning: Some nans were detected.\n'); end
		else
			Part_AUC = nan(n_study, n_rep);
			for si=1:n_study
				for ri=1:n_rep
					res_name = sprintf('%sDID_CVT01_Si%02d-Ri%03d_%s_*_MSN-%03d_MTN-%s.mat', result_path, si, ri, Net_Name, MAX_N_SN, Method_Name);
					file_info = dir(res_name);
					if numel(file_info)~=1
						fprintf('************ Ignoring [%60s]\n', res_name);
						continue;
					end
					fprintf('Reading from: [%s]\n', file_info.name);
					res_data = load([result_path file_info.name]);
					Part_AUC(si, ri) = res_data.te_auc;
				end
			end
			save(sav_name, 'Part_AUC', 'Method_Name', 'Net_Name', 'MAX_N_SN');
		end
		Result_AUC(mi, ni, :, :) = Part_AUC;
	end
end
if any(isnan(Result_AUC(:)))
	fprintf('[i] Warning some values are missing.\n');
end

%% Saving
save('./Collected_NOP_AUCs.mat', 'Result_AUC', 'net_lst', 'method_lst');
