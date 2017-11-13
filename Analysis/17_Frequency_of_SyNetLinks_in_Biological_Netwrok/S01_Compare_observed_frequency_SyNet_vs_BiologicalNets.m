clc;
clear;

%% Initialization
addpath('../_Utilities/');
run_id = 100;
n_lnk = 1000;
n_rep = 10000;
net_name = 'STRING';
net_opt.GE_Path = getPath('SyNet');
ge_data = load(net_opt.GE_Path, 'Gene_Name');
net_opt.PreferredGenes = ge_data.Gene_Name;
net_opt.MAX_N_PAIR = 10000;

%% Load network
fprintf('Loading [%s] network.\n', net_name);
net_info = LoadNetworkAdj(net_name, net_opt);
Net_Adj = net_info.Net_Adj;
Net_GeneName = net_info.Gene_Name;
clear net_info
n_gene = numel(Net_GeneName);
if min(Net_Adj(:))<0, error('Not implemented for negative links'); end

%% Load SyNet
SyNet_path = '../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat';
SyNet_info = load(SyNet_path, 'PP_Info', 'Gene_Name');
SyNet_lnk = SyNet_info.PP_Info(1:10000,:);
SyNet_GeneName = SyNet_info.Gene_Name;
SyNet_Map = containers.Map();
for si=1:size(SyNet_lnk,1)
    Pair_str = [SyNet_GeneName{SyNet_lnk(si,1)} ';' SyNet_GeneName{SyNet_lnk(si,2)}];
    SyNet_Map(Pair_str) = si;
end

%% Make reference gene set
InSyNet = ismember(Net_GeneName, SyNet_GeneName);

%% Main loop
fprintf('[%d] Random selection of [%d] links from [%s] ... \n', n_rep, n_lnk, net_name);
OL_Freq = zeros(n_rep, 1);
for ri=1:n_rep
    showprogress(ri, n_rep);
    
    %% Collect N random links
    Net_Lnk = [];
    while size(Net_Lnk,1)<n_lnk
        rnd_ind = randi(n_gene, n_lnk*10, 2);
        pair_ind = sub2ind([n_gene n_gene], rnd_ind(:,1), rnd_ind(:,2));
        has_lnk = Net_Adj(pair_ind)>0;
        Net_Lnk = [Net_Lnk; rnd_ind(has_lnk,:)];
    end
    Net_Lnk = Net_Lnk(1:n_lnk, :);
    
    %% Filter non-existing genes
    is_val = InSyNet(Net_Lnk);
    
    %% Measure overlap
    for li=1:n_lnk
        Pair_str1 = [Net_GeneName{Net_Lnk(li,1)} ';' Net_GeneName{Net_Lnk(li,2)}];
        Pair_str2 = [Net_GeneName{Net_Lnk(li,2)} ';' Net_GeneName{Net_Lnk(li,1)}];
        if SyNet_Map.isKey(Pair_str1) || SyNet_Map.isKey(Pair_str2)
            OL_Freq(ri) = OL_Freq(ri) + 1;
        end
    end
end

%% Save output
sav_name = sprinf('./Net_SyNet_Overlap/Net-OV_%s_NL%d_%d', net_name, n_lnk, run_id);
fprintf('Saving the results in [%s]\n', sav_name);
save(sav_name, 'OL_Freq', 'net_name', 'Gene_Name', 'SyNet_path');




