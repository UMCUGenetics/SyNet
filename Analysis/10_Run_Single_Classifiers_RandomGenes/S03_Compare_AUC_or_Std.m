clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/BoxplotEx');
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');

%% Loading results
Method_lst = {};
Net_lst = {};
auc_mat = [];
std_mat = [];
for cv_id=[1 91]
	res_name = sprintf('Collected_NOP_AUCs_CV%02d.mat', cv_id);
	fprintf('Loading results from [%s] ...\n', res_name);
	clc_data = load(res_name);
	Method_lst = [Method_lst strcat(num2str(cv_id), '-', clc_data.method_lst)];
	Net_lst = clc_data.net_lst;
	res_auc = squeeze(clc_data.Result_AUC);
	res_auc = median(res_auc, 4);
	std_mat = [std_mat; std(res_auc,0,3)];
	auc_mat = [auc_mat; mean(res_auc,3)];
end
n_met = numel(Method_lst);
n_net = numel(Net_lst);
n_std = 14;
n_rep = 10;

imagesc(auc_mat);
colormap(jet(10));
% colormap(copper(10));
set(gca, 'XTick', 1:n_net, 'XTickLabel', Net_lst, 'XTickLabelRotation', 20, ...
	'YTick', 1:n_met, 'YTickLabel', Method_lst);
colorbar();
% caxis([0.63 0.68]);

figure();
auc_mat = auc_mat([1], :);
std_mat = std_mat([1], :);
bar(auc_mat, 'Grouped');
colormap(jet(n_net));
ylim([0.55 0.65]);

