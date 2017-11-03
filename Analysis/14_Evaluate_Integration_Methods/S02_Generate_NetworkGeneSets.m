clc;
clear;

%% Initialization
addpath('../11_Perform_LassoTypes/');
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
addpath('../_Utilities/');
sav_path = 'NetNei_Files/';
[~,~] = mkdir(sav_path);
MAX_SUBNET_SIZE = 20;
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name');
n_gene = numel(Gene_Name);
GMap = containers.Map(Gene_Name, 1:n_gene);
net_name = 'STRING';
% net_name = 'Random';

%% Load network
switch net_name
    case 'Random'
        n_gene = numel(Gene_Name);
        Net_Adj = rand(n_gene);
    case 'STRING'
        net_path = getPath(net_name);
		fid = fopen(net_path, 'r');
		Header_lst = regexp(fgetl(fid), '\t', 'split');
		if numel(Header_lst)==2
			fprintf('No weight exists. Selecting all links.\n');
			net_cell = textscan(fid, '%s%s', 'Delimiter', '\t', 'ReturnOnError', 0);
			if ~feof(fid), error(); end
		else
			fprintf('Selecting of links from top weighted interactions.\n');
			net_cell = textscan(fid, '%s%s%d', 'Delimiter', '\t', 'ReturnOnError', 0);
            [vid, sid] = sort(net_cell{3}, 'Descend');
            net_cell = {net_cell{1}(sid) net_cell{2}(sid)};
		end
		fclose(fid);
		n_lnk = numel(net_cell{1});
        if n_lnk > realmax('single'), error(); end
		fprintf('[i] Network contains [%d] links.\n', n_lnk);
		Net_Adj = zeros(n_gene, 'single');
		fprintf('Forming the Adj matrix with [%d] genes: ', n_gene);
		for ii=1:n_lnk
			showprogress(ii, n_lnk);
			if GMap.isKey(net_cell{1}{ii}) && GMap.isKey(net_cell{2}{ii})
				gi = GMap(net_cell{1}{ii});
				gj = GMap(net_cell{2}{ii});
				Net_Adj(gi, gj) = n_lnk - ii + 1;
				Net_Adj(gj, gi) = n_lnk - ii + 1;
			end
		end
		clear net_cell
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
	SubNet_Full{gi} = getNetNeighborsBreadthFirst(Neig_cell, Neig_cell{gi}, MAX_SUBNET_SIZE+1, 1);
end

%% Save network neighbors
sav_name = sprintf([sav_path 'NetNei_%s_NN%02d.mat'], net_name, MAX_SUBNET_SIZE);
fprintf('Saving network file to [%s]\n', sav_name);
save(sav_name, 'SubNet_Full', 'Gene_Name');

