function S03_Overlap_SyNet_vs_Net_PerThresholded_Genes(Ref_Name, net_name, SHUFFLE_NODES)
%{
for ni in HumanInt BioPlex BioGRID IntAct STRING HBBrain HBKidney HBOvary HBLympNode HBGland; do
PARAM=\'SyNet\',\'$ni\',0;
echo $PARAM;
sbatch --job-name=GeneOL-$PARAM --output=Logs/GNOL-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=10GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S03_Overlap_SyNet_vs_Net_PerThresholded_Genes "$PARAM";
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
Ref_GeneName = SyNet_GeneName(unique(SyNet_lnk', 'Stable'));
clear SyNet_GeneName SyNet_lnk
fprintf('Reference network [%s] has [%d] genes.\n', Ref_Name, numel(Ref_GeneName));

%% Load network
fprintf('Loading [%s] network.\n', net_name);
net_info = LoadNetworkAdj(net_name, net_opt);
Net_Adj = net_info.Net_Adj;
Net_GeneName = net_info.Gene_Name;
if net_opt.Shuffle_Nodes
    net_name = [net_name '-SHFL'];
end
clear net_info
if min(Net_Adj(:))<0, error('Not implemented for negative links'); end
Net_nGene = numel(Net_GeneName);

%% Convert to pair indices
Net_Adj = triu(Net_Adj,1);
[Net_sval, Net_EdgeIndex] = sort(Net_Adj(:), 'Descend');
clear Net_Adj
Net_EdgeIndex = Net_EdgeIndex(Net_sval>0);
[Net_SelLnk(:,1), Net_SelLnk(:,2)] = ind2sub([Net_nGene Net_nGene], Net_EdgeIndex);
clear Net_sval Net_EdgeIndex
Top_Gene_lst = Net_GeneName(unique(Net_SelLnk', 'stable'));
if numel(Top_Gene_lst)~=Net_nGene, error(); end

%% Main loop
for ti=1:n_thresh
    n_SelGene = floor(Net_nGene * Ratio_lst(ti));
    fprintf('== [%04d/%04d] #Total genes=%d, selecting top [%d] genes ...\n', ti, n_thresh, Net_nGene, n_SelGene);
    
    Real_Freq = sum(ismember(Ref_GeneName, Top_Gene_lst(1:n_SelGene)));
    
    %% Randomized version
    Rand_Freq = zeros(n_rep, 1);
    for ri=1:n_rep
        if mod(ri,200)==0
            fprintf('\tRunning epoch [%5d/%5d] ...\n', ri, n_rep);
        end
        Rnd_GLst = Top_Gene_lst(randperm(Net_nGene));
        Rand_Freq(ri) = sum(ismember(Ref_GeneName, Rnd_GLst(1:n_SelGene)));
    end
    fprintf('Result: Read freq [%d], Mean Rnd [%0.2f/%0.2f] Std [%0.2f]\n', Real_Freq, mean(Rand_Freq), median(Rand_Freq), std(Rand_Freq));
    
    %% Save output
    sav_name = sprintf('./SyNet_OverlapPerThresh/GeneOV-PerThresh_%s_%s_MP%d_All_RV%0.2f.mat', Ref_Name, net_name, net_opt.MAX_N_PAIR, Ratio_lst(ti));
    fprintf('Saving the results in [%s]\n\n', sav_name);
    save(sav_name, 'Real_Freq', 'Rand_Freq', 'Ref_GeneName', 'Net_GeneName', 'Net_nGene', 'net_opt', 'n_SelGene', 'Ratio_lst');
end
end
