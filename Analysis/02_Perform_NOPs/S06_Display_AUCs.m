function S06_Display_AUCs()
clc;
clear;

%% Initialization
addpath('../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../Useful_Sample_Codes/getTop/');
addpath('../../Useful_Sample_Codes/ShowProgress/');
load('./Collected_Results.mat', 'res_info');
net_lst = {'Corr', 'Random', 'String', 'Multinet', 'KEGG', 'MSigDB', ...
	'DSN-HAIBE-T10000', 'DSN-META-T10000', 'DSN-ACES-T10000','DSN-HAIBERep-T10000', ...
	'CPN-HAIBE-T10000', 'CPN-META-T10000', 'CPN-ACES-T10000'};
% 	'DSN-HAIBERep-P99.90','DSN-HAIBERep-P99.99','DSN-HAIBERep-P99.999', ...
% 	'DSN-HAIBERep-T05000','DSN-HAIBERep-T10000','DSN-HAIBERep-T20000', ...
n_net = numel(net_lst);
in = ismember({res_info.Net_Name}, net_lst);
res_info = res_info(in);

%% Get data frequency
[comp_lst, comp_freq] = getTop({res_info.Data_Name});

%% Collect AUCs
cv_lst = unique({res_info.CV_ID});
n_cv = numel(cv_lst);
met_lst = unique({res_info.Method_Name});
n_met = numel(met_lst);
auc_mat = nan(n_net, n_met, n_cv);
all_cv = {res_info.CV_ID};
all_net = {res_info.Net_Name};
all_met = {res_info.Method_Name};
for ci=1:n_cv
	for ni=1:n_net
		for mi=1:n_met
			in_cv = ismember(all_cv, cv_lst{ci});
			in_net = ismember(all_net, net_lst{ni});
			in_met = ismember(all_met, met_lst{mi});
			in = in_cv & in_net & in_met;
			if any(in)
				if ~isnan(auc_mat(ni, mi, ci))
					error();
				else
					auc_mat(ni, mi, ci) = res_info(in).TE_AUC;
				end
			end
		end
	end
end
avg_mat = mean(auc_mat, 3, 'omitnan');
avg_mat(end+1, :) = mean(avg_mat, 1, 'omitnan');
avg_mat(:, end+1) = mean(avg_mat, 2, 'omitnan');
disp(array2table(avg_mat, 'RowNames', [net_lst 'AVG'] , 'VariableNames', [met_lst 'AVG']));
fprintf([repmat('$',1,50) '\n\n']);

%% Improvement per Methods
close all
avg_net = zeros(n_net, n_met);
std_net = zeros(n_net, n_met);
for mi=1:n_met
	for ni=1:n_net
		avg_net(ni,mi) = mean(auc_mat(ni,mi,:), 3, 'omitnan');
		std_net(ni,mi) = std(auc_mat(ni,mi,:), 0, 3, 'omitnan');
	end
end
figure('Position', [100 100 1200 400]);
bar(avg_net', 'grouped');
legend(net_lst, 'Location', 'Best');
set(gca, 'XTickLabel', met_lst);
ylim([0.62 0.75]);
colormap(lines);

%% Saving
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [12 4], 'PaperPosition', [0 0 12 4]);
print('-dpdf', '-r300', './Plots/S06_OverallPerformance_PerMethod.pdf');
end

