clc;
clear;

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
addpath('../19_Compute_Topological_Measures_For_Network/');
Net_lst = {'STRING' 'KEGG' 'MSigDB' 'HPRD' 'I2D' 'IntAct' 'HumanInt' 'BioPlex' 'BioGRID'};
n_net = numel(Net_lst);
net_opt.MAX_N_PAIR = 50000;

%% Loading networks and unify gene list
RefGene_lst = {};
for ni=1:n_net
    fprintf('Loading [%s] network\n', Net_lst{ni});
    net_info = LoadNetworkAdj(Net_lst{ni}, net_opt);
    if ni==1
        RefGene_lst = net_info.Gene_Name;
    else
        RefGene_lst = intersect(RefGene_lst, net_info.Gene_Name);
    end
end

%% Save the unified gene list

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
