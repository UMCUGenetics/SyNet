function S01_EvaluateNetOverlap(Ref_Name, net_name, SHUFFLE)
%% Run: PARAM="'AvgSynACr','STRING',0"; sbatch --job-name=NO-$PARAM --output=Logs/NO-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=5GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S01_EvaluateNetOverlap "$PARAM";
clc;

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath('../_Utilities/');
net_opt.GE_Path = getPath('SyNet');
ge_data = load(net_opt.GE_Path, 'Gene_Name');
net_opt.PreferredGenes = ge_data.Gene_Name;
net_opt.MAX_N_PAIR = 25000;
SampleSize = 10000;
n_rep = 10000;
if ismac
    Ref_Name = 'SyNet';
    net_name = 'HBBlood';
    %net_name = 'AbsCorr';
    SHUFFLE = 0;
end

%% Load SyNet
if strcmp(Ref_Name, 'SyNet')
    SyNet_path = '../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat';
    SyNet_info = load(SyNet_path, 'PP_Info', 'Gene_Name');
    SyNet_lnk = SyNet_info.PP_Info(1:10000,1:2);
    SyNet_GeneName = SyNet_info.Gene_Name;
else
    net_info = LoadNetworkAdj(Ref_Name);
    SyNet_GeneName = net_info.Gene_Name;
    n_gene = numel(SyNet_GeneName);
    Tmp_Adj = triu(net_info.Net_Adj, 1);
    [~, s_ind] = sort(Tmp_Adj(:), 'Descend');
    [SyNet_lnk(:,1), SyNet_lnk(:,2)] = ind2sub([n_gene n_gene], s_ind(1:10000));
    clear net_info Tmp_Adj s_ind
end
SyNet_Map = [
    strcat(SyNet_GeneName(SyNet_lnk(:,1)), ';', SyNet_GeneName(SyNet_lnk(:,2)));
    strcat(SyNet_GeneName(SyNet_lnk(:,2)), ';', SyNet_GeneName(SyNet_lnk(:,1)))
    ];

%% Load network
fprintf('Loading [%s] network.\n', net_name);
net_info = LoadNetworkAdj(net_name, net_opt);
Net_Adj = net_info.Net_Adj;
Net_GeneName = net_info.Gene_Name;
clear net_info
if min(Net_Adj(:))<0, error('Not implemented for negative links'); end
Net_Adj = triu(Net_Adj,1);
Net_EdgeIndex = find(Net_Adj(:)>0);
clear Net_Adj
Net_nlnk = numel(Net_EdgeIndex);
Net_ngene = numel(Net_GeneName);
if SHUFFLE == 1
    net_name = [net_name '-SHFL'];
end

%% Main loop
fprintf('[%d] Random selection of [%d] links from [%s] ... \n', n_rep, SampleSize, net_name);
OL_Freq = zeros(n_rep, 1);
Net_SelLnk = zeros(SampleSize, 2);
for ri=1:n_rep
    if mod(ri,200)==0
        fprintf('Running epoch [%5d/%5d] ...\n', ri, n_rep);
    end
    
    %% Collect N random links
    Rnd_Index = randi(Net_nlnk, SampleSize, 1);
    [Net_SelLnk(:,1), Net_SelLnk(:,2)] = ind2sub([Net_ngene Net_ngene], Net_EdgeIndex(Rnd_Index));
    
    %% Measure overlap
    if SHUFFLE == 1
        Net_GeneName = Net_GeneName(randperm(Net_ngene));
    end
    Pair_str = strcat(Net_GeneName(Net_SelLnk(:,1)), ';', Net_GeneName(Net_SelLnk(:,2)));
    is_in = ismember(Pair_str, SyNet_Map);
    OL_Freq(ri) = sum(is_in);
end

%% Save output
sav_name = sprintf('./SyNet_Overlap/NetOV_%s_%s_NL%d.mat', Ref_Name, net_name, SampleSize);
fprintf('Saving the results in [%s]\n', sav_name);
save(sav_name, 'OL_Freq', 'net_name', 'Net_GeneName', 'SyNet_GeneName', 'Net_nlnk', 'net_opt', 'SampleSize');
end


