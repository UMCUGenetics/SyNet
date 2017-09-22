clc;
clear;

%% Initialization
result_path = '../08_Perform_NOPs_Restricted_Network/Results_Files/';
method_lst = {'TAgNMC','TNMC','TLEx','TAgLEx'};
n_met = numel(method_lst);
net_lst = {
	'Random-T00010', 'Random-T00020', 'Random-T00050', ...
	'Random-T00100', 'Random-T00200', 'Random-T00500', ...
	};
n_net = numel(net_lst);
cv_ind = 91;
n_study = 1;
n_rep = 50;
sav_path = './Collected_AUCs/';
[~,~] = mkdir(sav_path);

%% Main loop
Result_AUC = nan(n_met, n_net, n_study, n_rep);
for mi=1:n_met
	for ni=1:n_net
		Method_Name = method_lst{mi};
		Net_Name = net_lst{ni};
		sav_name = sprintf([sav_path 'CRES_CVT%02d_%s_%s.mat'], cv_ind, Method_Name, Net_Name);
		if exist(sav_name, 'file')
			fprintf('/// Loading from [%s]\n', sav_name);
			res_data = load(sav_name);
			if ~isequal(Method_Name, res_data.Method_Name) || ~isequal(Net_Name, res_data.Net_Name), error(); end
			Part_AUC = res_data.Part_AUC;
			if any(isnan(Part_AUC(:))), fprintf('Warning: Some nans were detected.\n'); end
		else
			Part_AUC = nan(n_study, n_rep);
			for si=1:n_study
				for ri=1:n_rep
					res_name = sprintf('%sDID_CVT%02d_Si%02d-Ri%03d_%s_*_MTN-%s.mat', result_path, cv_ind, si, ri, Net_Name, Method_Name);
					file_info = dir(res_name);
					if numel(file_info)~=1
						fprintf('************ Ignoring [%60s]\n', res_name);
						%error();
						continue;
					end
					fprintf('Reading from: [%s]\n', file_info.name);
					res_data = load([result_path file_info.name]);
					Part_AUC(si, ri) = res_data.te_auc;
				end
			end
			save(sav_name, 'Part_AUC', 'Method_Name', 'Net_Name');
		end
		Result_AUC(mi, ni, :, :) = Part_AUC;
	end
end
if any(isnan(Result_AUC(:)))
	fprintf('[i] Warning some values are missing.\n');
end

%% Saving
out_name = sprintf('./Collected_NOP_AUCs_CV%02d.mat', cv_ind);
save(out_name, 'Result_AUC', 'net_lst', 'method_lst');
