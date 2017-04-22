
%% Initialization
clc;
clear;
n_study = 14;
targ_ind = 13;
rest_ind = setdiff(1:n_study, targ_ind);

%% loading pairwise data
opw_data = load('./PWR_Files_CV01/SyNet/PWR_SyNet_00000001-00100000.mat');
old_auc = [opw_data.auc_pair(:,1:2) median(opw_data.auc_pair(:,rest_ind+2),2)];

npw_data = load('./PWR_Files/SyNet/PWR_SyNet_00000001-00030000.mat');
new_auc = [npw_data.auc_pair(:,1:2) median(npw_data.auc_pair(:,rest_ind+2),2)];
npw_data = load('./PWR_Files/SyNet/PWR_SyNet_00030001-00060000.mat');
new_auc = [new_auc; npw_data.auc_pair(:,1:2) median(npw_data.auc_pair(:,rest_ind+2),2)];
npw_data = load('./PWR_Files/SyNet/PWR_SyNet_00060001-00090000.mat');
new_auc = [new_auc; npw_data.auc_pair(:,1:2) median(npw_data.auc_pair(:,rest_ind+2),2)];
n_pair = size(new_auc, 1);
old_auc = old_auc(1:n_pair,:);
if ~isequal(old_auc(:,1:2), new_auc(:,1:2))
	error();
end

%% Calculate performance
addpath('../../../../Useful_Sample_Codes/getAUC/');
addpath('../02_Perform_NOPs');
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Expression', 'Patient_Label', 'Study_Index');
zData = zscore(Gene_Expression);
auc_mat = [old_auc new_auc(:,3) zeros(n_pair,1)];
in = Study_Index==targ_ind;
for pi=1:n_pair
	showprogress(pi, n_pair, 50);
	mData = mean(zData(in, auc_mat(pi,1:2)),2);
	auc_mat(pi,5) = measureAUC(mData, Patient_Label(in), 20);
end

%% Sorting
[~, sid] = sort(auc_mat(:,5), 'Descend');
auc_mat = auc_mat(sid,:);

%% Plotting
% hold on
% plot(auc_mat(:,3), 'k');
% plot(auc_mat(:,4), 'b');
% plot(auc_mat(:,5), 'r');
corr(auc_mat(:,3:5), 'Type', 'Spearman')
% legend({'Old' 'New' 'Test'});
