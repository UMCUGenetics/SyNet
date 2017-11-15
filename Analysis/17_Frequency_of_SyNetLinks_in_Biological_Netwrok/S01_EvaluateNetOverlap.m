function S01_EvaluateNetOverlap(net_name, SampleSize, n_rep, SHUFFLE)
%% Run: PARAM="'STRING',1000,10000,0"; sbatch --job-name=NO-$PARAM --output=Logs/NO-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=7GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S01_EvaluateNetOverlap "$PARAM";

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath('../_Utilities/');
net_opt.GE_Path = getPath('SyNet');
ge_data = load(net_opt.GE_Path, 'Gene_Name');
net_opt.PreferredGenes = ge_data.Gene_Name;
net_opt.MAX_N_PAIR = 50000;
if ismac
    SampleSize = 1000;
    n_rep = 10000;
    net_name = 'STRING';
    net_name = 'AbsCorr';
    SHUFFLE = 1;
end

%% Load network
fprintf('Loading [%s] network.\n', net_name);
net_info = LoadNetworkAdj(net_name, net_opt);
Net_Adj = net_info.Net_Adj;
Net_GeneName = net_info.Gene_Name;
clear net_info
n_gene = numel(Net_GeneName);
if min(Net_Adj(:))<0, error('Not implemented for negative links'); end
if SHUFFLE == 1
    net_name = [net_name '-SHFL'];
end
Net_nlnk = numel(nonzeros(triu(Net_Adj)));

%% Load SyNet
SyNet_path = '../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet.mat';
SyNet_info = load(SyNet_path, 'PP_Info', 'Gene_Name');
SyNet_lnk = SyNet_info.PP_Info(1:10000,:);
SyNet_GeneName = SyNet_info.Gene_Name;
SyNet_Map = containers.Map();
for si=1:size(SyNet_lnk,1)
    Pair_str = [SyNet_GeneName{SyNet_lnk(si,1)} ';' SyNet_GeneName{SyNet_lnk(si,2)}];
    SyNet_Map(Pair_str) = si;
    Pair_str = [SyNet_GeneName{SyNet_lnk(si,2)} ';' SyNet_GeneName{SyNet_lnk(si,1)}];
    SyNet_Map(Pair_str) = si;
end

%% Main loop
fprintf('[%d] Random selection of [%d] links from [%s] ... \n', n_rep, SampleSize, net_name);
OL_Freq = zeros(n_rep, 1);
for ri=1:n_rep
    showprogress(ri, n_rep);
    
    %% Collect N random links
    Net_Lnk = [];
    while size(Net_Lnk,1)<SampleSize
        rnd_ind = randi(n_gene, SampleSize*50, 2);
        pair_ind = sub2ind([n_gene n_gene], rnd_ind(:,1), rnd_ind(:,2));
        has_lnk = Net_Adj(pair_ind)>0;
        Net_Lnk = [Net_Lnk; rnd_ind(has_lnk,:)];
    end
    Net_Lnk = Net_Lnk(1:SampleSize, :);
    
    if SHUFFLE == 1
        Net_GeneName = Net_GeneName(randperm(n_gene));
    end
    
    %% Measure overlap
    for li=1:SampleSize
        Pair_str = [Net_GeneName{Net_Lnk(li,1)} ';' Net_GeneName{Net_Lnk(li,2)}];
        if SyNet_Map.isKey(Pair_str)
            OL_Freq(ri) = OL_Freq(ri) + 1;
        end
    end
end

%% Plotting
% plot(sort(OL_Freq));

%% Save output
sav_name = sprintf('./SyNet_Overlap/NetOV_%s_NL%d.mat', net_name, SampleSize);
fprintf('Saving the results in [%s]\n', sav_name);
save(sav_name, 'OL_Freq', 'net_name', 'Net_GeneName', 'SyNet_GeneName', 'SyNet_path', 'Net_nlnk', 'net_opt', 'SampleSize');
end


