clc;
clear;

%% Initialization
addpath('../11_Perform_LassoTypes/');
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
sav_path = 'NetNei_Files/';
[~,~] = mkdir(sav_path);
n_epoch = 10000;
MAX_SUBNET_SIZE = 5;
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name');
n_gene = numel(Gene_Name);
net_name = 'Random';

%% Load network
switch net_name
    case 'Random'
        Net_Adj = rand(n_gene);
    otherwise
        fprintf('Unknown method');
end

%% Generate Neighbor Sets
fprintf('Generating neighbor sets and subnetworks: \n');
Neig_cell = getNeighborsFromAdj(Net_Adj, MAX_SUBNET_SIZE);

%% Generate Subnetworks
SubNet_Full = cell(n_gene, 1);
for gi=1:n_gene
	showprogress(gi, n_gene);
	SubNet_Full{gi} = getNetNeighborsBreadthFirst(Neig_cell, Neig_cell{gi}, MAX_SUBNET_SIZE, 1);
end

%% Save network neighbors
sav_name = sprintf([sav_path 'NetNei_%s_NN%02d.mat'], net_name, MAX_SUBNET_SIZE);
fprintf('Saving network file to [%s]\n', sav_name);
save(sav_name, 'SubNet_Full', 'Gene_Name');

