function te_auc = S04_Display_CVAUC(cv_id)
clc;
% clear;

%% Initialization
result_path = './Results_Files/';
net_lst = {'DSN-ACES-T05000', 'DSN-ACES-T10000', 'DSN-ACES-T20000', 'DSN-ACES-P99.99', 'DSN-ACES-P99.9', 'DSN-ACES-P99', 'String', 'Corr', 'Random'};
n_net = numel(net_lst);
method_lst = {'iPark', 'iChuang', 'iTaylor', 'Feral'};
n_met = numel(method_lst);

te_auc = nan(n_net, n_met);
for ni=1:n_net
	for mi=1:n_met
		result_temp = sprintf([result_path 'CID-' cv_id '*_NET-' net_lst{ni} '_*MSN-100_MTN-' method_lst{mi} '.mat']);
		ds_info = dir(result_temp);
		if numel(ds_info)==0
			fprintf('Ignoring [%s]\n', result_temp);
			continue;
		end
		dataset_name = ds_info.name;
		fprintf('Reading from: [%s]\n', ds_info.name);
		data = load([result_path ds_info(1).name], 'te_auc');
		te_auc(ni, mi) = round(data.te_auc*100, 3);
	end
end
te_auc(end+1, :) = mean(te_auc,1, 'omitnan');
te_auc(:, end+1) = mean(te_auc,2, 'omitnan');

fprintf('\n\n\nDataset: %s\n', dataset_name);
disp(array2table(te_auc, 'RowNames', [net_lst 'AVG'] , 'VariableNames', [method_lst 'AVG']));
end