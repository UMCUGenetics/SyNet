function S04_Overlap_SyNet_vs_Net_PerThresholded_Links(Ref_Name, net_name, LIMIT_GENES, SHUFFLE_NODES)
%{
for mi in 0 1; do
for ni in HumanInt BioPlex BioGRID IntAct STRING HBBrain HBKidney HBOvary HBLympNode HBGland; do
PARAM=\'SyNet\',\'$ni\',$mi,0;
echo $PARAM;
sbatch --job-name=LnkOL-$PARAM --output=Logs/LnkOL-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=10GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S04_Overlap_SyNet_vs_Net_PerThresholded_Links "$PARAM";
done; 
read -p "`date`: $PARAM. Press a key" -t 180
done
%}
clc;

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath('../_Utilities/');
if ~exist('SHUFFLE_NODES', 'var'), SHUFFLE_NODES = 0; end
net_opt.GE_Path = getPath('SyNet');
ge_data = load(net_opt.GE_Path, 'Gene_Name');
net_opt.PreferredGenes = ge_data.Gene_Name;
net_opt.MAX_N_PAIR = 50000;
net_opt.Shuffle_Nodes = SHUFFLE_NODES;
N_Ref_lnk = 3544;
n_rep = 1000;
Ratio_lst = 0.05:0.05:1;
n_thresh = numel(Ratio_lst);
if ismac
    Ref_Name = 'SyNet'; %'AvgSyn'
    net_name = 'HumanInt';
    LIMIT_GENES = 0;
end

%% Load SyNet
if strcmp(Ref_Name, 'SyNet')
    SyNet_path = '../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet_AvgSynACr.mat';
    SyNet_info = load(SyNet_path, 'PP_Info', 'Gene_Name');
    SyNet_lnk = SyNet_info.PP_Info(1:N_Ref_lnk,1:2);
    SyNet_GeneName = SyNet_info.Gene_Name;
    clear SyNet_info
else
    net_info = LoadNetworkAdj(Ref_Name);
    SyNet_GeneName = net_info.Gene_Name;
    n_gene = numel(SyNet_GeneName);
    Tmp_Adj = triu(net_info.Net_Adj, 1);
    [~, s_ind] = sort(Tmp_Adj(:), 'Descend');
    [SyNet_lnk(:,1), SyNet_lnk(:,2)] = ind2sub([n_gene n_gene], s_ind(1:N_Ref_lnk));
    clear net_info Tmp_Adj s_ind
end
SyNet_Map = [
    strcat(SyNet_GeneName(SyNet_lnk(:,1)), ';', SyNet_GeneName(SyNet_lnk(:,2)));
    strcat(SyNet_GeneName(SyNet_lnk(:,2)), ';', SyNet_GeneName(SyNet_lnk(:,1)))
    ];
fprintf('Reference network [%s] has [%d] genes and [%d] pairs.\n', Ref_Name, numel(unique(SyNet_lnk)), size(SyNet_lnk,1));

%% Load network
fprintf('Loading [%s] network.\n', net_name);
if LIMIT_GENES
    net_opt.PreferredGenes = SyNet_GeneName(unique(SyNet_lnk(:)));
    limit_method = 'LimitedToRef';
else
    limit_method = 'All';
end
net_info = LoadNetworkAdj(net_name, net_opt);
Net_Adj = net_info.Net_Adj;
Net_GeneName = net_info.Gene_Name;
if net_opt.Shuffle_Nodes
    net_name = [net_name '-SHFL'];
end
clear net_info
if min(Net_Adj(:))<0, error('Not implemented for negative links'); end
Net_ngene = numel(Net_GeneName);

%% Convert to pair indices
Net_Adj = triu(Net_Adj,1);
[Net_sval, Net_EdgeIndex] = sort(Net_Adj(:), 'Descend');
clear Net_Adj
Net_EdgeIndex = Net_EdgeIndex(Net_sval>0);
clear Net_sval
Net_nlnk = numel(Net_EdgeIndex);

%% Main loop
for ti=1:n_thresh
    n_TopPair = floor(Net_nlnk * Ratio_lst(ti));
    fprintf('== [%04d/%04d] #Total Links=%d, selecting top [%d] links ...\n', ti, n_thresh, Net_nlnk, n_TopPair);
    
    Net_SelLnk = zeros(n_TopPair, 2);
    [Net_SelLnk(:,1), Net_SelLnk(:,2)] = ind2sub([Net_ngene Net_ngene], Net_EdgeIndex(1:n_TopPair));
    Real_Freq = getLinkOverlap(Net_SelLnk, Net_GeneName, SyNet_Map);
    
    %% Randomized version
    Rand_Freq = zeros(n_rep, 1);
    for ri=1:n_rep
        if mod(ri,200)==0
            fprintf('\tRunning epoch [%5d/%5d] ...\n', ri, n_rep);
        end
        Rand_Freq(ri) = getLinkOverlap(Net_SelLnk, Net_GeneName(randperm(Net_ngene)), SyNet_Map);
    end
    fprintf('Result: Read freq [%d], Mean Rnd [%0.2f/%0.2f] Std [%0.2f]\n', Real_Freq, mean(Rand_Freq), median(Rand_Freq), std(Rand_Freq));
    
    %% Save output
    sav_name = sprintf('./SyNet_OverlapPerThresh/LnkOV-PerThresh_%s_%s_MP%d_%s_RV%0.2f.mat', Ref_Name, net_name, net_opt.MAX_N_PAIR, limit_method, Ratio_lst(ti));
    fprintf('Saving the results in [%s]\n\n', sav_name);
    save(sav_name, 'Real_Freq', 'Rand_Freq', 'Net_GeneName', 'SyNet_GeneName', 'Net_nlnk', 'net_opt', 'n_TopPair', 'Ratio_lst');
end
end

function overlap_freq = getLinkOverlap(Net_SelLnk, Net_GeneName, SyNet_Map)
Pair_str = strcat(Net_GeneName(Net_SelLnk(:,1)), ';', Net_GeneName(Net_SelLnk(:,2)));
is_in = ismember(Pair_str, SyNet_Map);
overlap_freq = sum(is_in);
end