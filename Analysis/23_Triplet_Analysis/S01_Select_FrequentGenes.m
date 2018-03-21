clc;
clear;

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
addpath('../19_Compute_Topological_Measures_For_Network/');
Net_lst = {'STRING' 'KEGG' 'MSigDB' 'HPRD' 'I2D' 'IntAct' 'HumanInt' 'BioPlex' 'BioGRID'};
n_net = numel(Net_lst);
MAX_N_PAIR = 50000;

%% Load GE data
GeneExpression_Path = getPath('SyNet');
GE_Info = load(GeneExpression_Path, 'Gene_Name');

%% Loading networks and unify gene list
GMap = containers.Map();
Gene_InNet = zeros(2e4, n_net);
Gene_Name = cell(2e4, 1);
net_opt.PreferredGenes = GE_Info.Gene_Name;
net_opt.MAX_N_PAIR = MAX_N_PAIR;
for ni=1:n_net
    fprintf('Loading [%s] network\n', Net_lst{ni});
    net_info = LoadNetworkAdj(Net_lst{ni}, net_opt);
    Gene_Lst = net_info.Gene_Name;
    for gi=1:numel(Gene_Lst)
        if GMap.isKey(Gene_Lst{gi})
            g_ind = GMap(Gene_Lst{gi});
        else
            g_ind = GMap.Count + 1;
            GMap(Gene_Lst{gi}) = g_ind;
            Gene_Name{g_ind} = Gene_Lst{gi};
        end
        Gene_InNet(g_ind, ni) = 1;
    end
end
Gene_InNet(GMap.Count+1:end, :) = [];
Gene_Name(GMap.Count+1:end) = [];
clear Gene_Lst

%% Sorting gene set
Gene_InFreq = sum(Gene_InNet, 2);
[~, sid] = sort(Gene_InFreq, 'Descend');
Gene_InNet = Gene_InNet(sid, :);
Gene_Name = Gene_Name(sid);

%% Plotting
imagesc(Gene_InNet);
set(gca, 'XTick', 1:n_net, 'XTickLabel', Net_lst);

%% Select top genes
Ref_GeneName = Gene_Name(Gene_InFreq>=n_net-1);
n_ref = numel(Ref_GeneName);
Ref_GeneIndex = zeros(n_ref, 1);
for ri=1:n_ref
    Ref_GeneIndex(ri, 1) = find(strcmp(GE_Info.Gene_Name, Ref_GeneName{ri}));
end

%% Save the unified gene list
save('./Gene_List/Reference_GList.mat', 'Ref_GeneName', 'Ref_GeneIndex', 'Gene_Name', 'Gene_InNet', 'Net_lst', 'MAX_N_PAIR');


