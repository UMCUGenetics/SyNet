function te_auc = S04_Collect_Test_AUC(cv_id)
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

if ~ismac
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
					data = load([result_path ds_info.name], 'te_auc');
					te_auc(mi,ni,si,ri) = round(data.te_auc*100, 3);
				end
			end
		end
	end
	save('./Collected_result.mat', 'te_auc');
else
	load('./Collected_result.mat');
	if any(isnan(te_auc(:))), error(); end
	load('/Users/amin/Technical/My_Works/Deep_Learning/113_Organize_SPADE_Codes/Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name', 'Study_Index');
	
	tmp = median(te_auc,4);
	tmp = squeeze(median(tmp, 3));
	boxplot(tmp);
	set(gca, 'XTick', 1:n_net, 'XTickLabel', net_lst, 'XTickLabelRotation', 25);
	
	auc_feral = squeeze(tmp(4,:,:));
	auc_single = squeeze(tmp(5,:,:));
	imagesc(auc_feral-auc_single)
	caxis([-6 6]); 
	colormap(flipud(jet(10)));
	colorbar();
	set(gca, 'YTick', 1:n_net, 'YTickLabel', net_lst, 'YTickLabelRotation', 25, ...
			 'XTick', 1:14, 'XTickLabel', Study_Name, 'XTickLabelRotation', 25);
	
	tmp = arrayfun(@(i) squeeze(te_auc(i,1,13,:)), 1:5, 'UniformOutput', false);
	tmp = [tmp{:}];
	
	close all
	figure('Position', [100 100 1200 500]);
	plot(tmp', 'LineWidth', 2);
	legend(method_lst);
	set(gca, 'XTick', 1:n_net, 'XTickLabel', net_lst, 'XTickLabelRotation', 25);
	
	boxplot(squeeze(median(median(te_auc,4),3)))
	
	tmp = squeeze(median(te_auc,4));
	tmp = squeeze(median(tmp,1));
	boxplot(tmp);
	set(gca, 'XTick', 1:n_net, 'XTickLabel', net_lst, 'XTickLabelRotation', 25);
end

% te_auc(end+1, :) = mean(te_auc,1, 'omitnan');
% te_auc(:, end+1) = mean(te_auc,2, 'omitnan');
%
% fprintf('\n\n\nDataset: %s\n', dataset_name);
% disp(array2table(te_auc, 'RowNames', [net_lst 'AVG'] , 'VariableNames', [method_lst 'AVG']));
end