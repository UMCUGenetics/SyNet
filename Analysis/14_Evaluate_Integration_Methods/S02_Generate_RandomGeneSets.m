clc;
clear;

%% Initialization
addpath('../11_Perform_LassoTypes/');
sav_path = 'NetNei_Files/';
[~,~] = mkdir(sav_path);
n_epoch = 10000;
MAX_SUBNET_SIZE = 20;
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name');
n_gene = numel(Gene_Name);
net_name = 'Random';

%% Generate Subnetworks
SubNet_Full = cell(n_epoch, 1);
for ei=1:n_epoch
	SubNet_Full{ei} = randperm(n_gene, MAX_SUBNET_SIZE);
end

%% Save network neighbors
sav_name = sprintf([sav_path 'NetNei_%s-NN%02d.mat'], net_name, MAX_SUBNET_SIZE);
fprintf('Saving network file to [%s]\n', sav_name);
save(sav_name, 'SubNet_Full', 'Gene_Name');

