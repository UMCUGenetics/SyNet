function S01_EvaluateNetOverlap(Ref_Name, net_name, SHUFFLE, LIMIT_GENES, MAX_N_PAIR)
%{ 
for ni in HumanInt BioPlex BioGRID IntAct STRING HBBrain HBKidney HBOvary HBLympNode HBGland; do 
for mi in 0 1; do  
PARAM=\'SyNet\',\'$ni\',$mi,0,100000; 
echo $PARAM; 
sbatch --job-name=NO-$PARAM --output=Logs/NO-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=10GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S01_EvaluateNetOverlap "$PARAM"; 
done; done
%}
clc;

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath('../_Utilities/');
if ismac
    Ref_Name = 'SyNet'; %'AvgSyn'
    net_name = 'STRING';
    %net_name = 'AbsCorr';
    SHUFFLE = 0;
    LIMIT_GENES = 1;
    MAX_N_PAIR = 50000;
end

net_opt.GE_Path = getPath('SyNet');
ge_data = load(net_opt.GE_Path, 'Gene_Name');
net_opt.PreferredGenes = ge_data.Gene_Name;
net_opt.MAX_N_PAIR = MAX_N_PAIR;
SampleSize = 3544;
N_Ref_lnk = 3544;
n_rep = 1000;

%% Load SyNet
if strcmp(Ref_Name, 'SyNet')
    SyNet_path = '../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet_AvgSynACr.mat';
    SyNet_info = load(SyNet_path, 'PP_Info', 'Gene_Name');
    SyNet_lnk = SyNet_info.PP_Info(1:N_Ref_lnk,1:2);
    SyNet_GeneName = SyNet_info.Gene_Name;
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
sav_name = sprintf('./SyNet_Overlap/NetOV_%s_%s_MP%d_%s_SS%d.mat', Ref_Name, net_name, net_opt.MAX_N_PAIR, limit_method, SampleSize);
fprintf('Saving the results in [%s]\n', sav_name);
save(sav_name, 'OL_Freq', 'net_name', 'Net_GeneName', 'SyNet_GeneName', 'Net_nlnk', 'net_opt', 'SampleSize');
end


