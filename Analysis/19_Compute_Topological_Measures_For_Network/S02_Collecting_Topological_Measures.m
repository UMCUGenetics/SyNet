clc;
clear;

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
% MAX_SyNet_Pair = 3544;
MAX_SyNet_Pair = 100000;
n_lnk = 50000;
Ref_Name = 'AvgSynACr';
% Shuff_Method = 'LnkShuff';
Shuff_Method = 'OneGRND';

%% Load labels
PI_Name = sprintf('./Topological_Data/PairInfo-%s_%s_MP%06d.mat', Shuff_Method, Ref_Name, MAX_SyNet_Pair);
fprintf('Loading pair info file [%s]\n', PI_Name);
Pair_Data = load(PI_Name);
Ref_GeneName = Pair_Data.SyNet_info.Gene_Name;
n_RefGeneName = numel(Ref_GeneName);
Pair_Info = [
    Pair_Data.PP_Info(:,1:2)  ones(MAX_SyNet_Pair, 1);
    Pair_Data.NP_Info(:,1:2) -ones(MAX_SyNet_Pair, 1);
    ];
TM_Label = double(Pair_Info(:,3));
n_Sample = size(Pair_Info, 1);
% clear Pair_Data

%% Collect data
net_lst = {'HumanInt' 'BioPlex','BioGRID','IntAct','STRING','HBBrain','HBKidney','HBOvary','HBLympNode','HBGland'}; % No 'AbsCorr'! as its used in SyNet
n_net = numel(net_lst);
tm_lst = {'DirectConnection' 'Jaccard' 'ClusteringCoefficient' 'ShortestPath' 'PageRank-FB0.95' 'PageRank-FB0.85' 'PageRank-FB0.75' 'PageRank-FB0.65' 'Eigenvector' 'Degree' 'Closeness' 'Betweenness'};
n_tm = numel(tm_lst);
TM_Data = zeros(MAX_SyNet_Pair*2, n_net*n_tm*2);
TM_Name = cell(n_net*n_tm*2, 1);
step = 1;
for ti=1:n_tm
    for ni=1:n_net
        data_name = sprintf('./Topological_Results/TM-%s_%s_%s_NP%06d_MSP%06d_%s.mat', Shuff_Method, Ref_Name, net_lst{ni}, n_lnk, MAX_SyNet_Pair, tm_lst{ti});
        fprintf('Loading feature from [%s].\n', data_name);
        Data_Info = load(data_name);
        if ~isequal(Pair_Data.SyNet_info.Gene_Name, Data_Info.Ref_GeneName), error(); end
        if ~isequal(Pair_Info(:,1:2), Data_Info.Pair_Info), error(); end
        
        TM_Name{step}   = sprintf('%s-%s-A', net_lst{ni}, tm_lst{ti});
        TM_Data(:,step)   = Data_Info.Pair_AvgScore;
        step = step + 1;
        
        if ismember(tm_lst{ti}, {'DirectConnection' 'ShortestPath' 'Jaccard'})
            if ~isequal(TM_Data(:,step-1), Data_Info.Pair_DifScore)
                error();
            end
        else
            TM_Name{step} = sprintf('%s-%s-D', net_lst{ni}, tm_lst{ti});
            TM_Data(:,step) = Data_Info.Pair_DifScore;
            step = step + 1;
        end
    end
end
TM_Data(:, step:end) = [];
TM_Name(step:end) = [];
[n_sample, n_feature] = size(TM_Data);

%% Normalization
%imagesc(TM_Data);
TM_Data_filtered = TM_Data;
for fi=1:n_feature
    is_inf = TM_Data_filtered(:,fi)==inf;
    if sum(is_inf)>n_sample*0.50
        fprintf('[w] Warning: Large number of samples in [%s] are inf [%d/%d] (%0.2f%%).\n', TM_Name{fi}, sum(is_inf), n_sample, sum(is_inf)*100/n_sample);
    end
    TM_Data_filtered(is_inf,fi) = max(TM_Data_filtered(~is_inf,fi))*1.01;
end
if any(isnan(TM_Data_filtered(:))) || any(TM_Data_filtered(:)==-inf), error(); end
TM_Data_z = zscore(TM_Data_filtered);
%imagesc(TM_Data_z);

%% Saving Data
sav_name = sprintf('./Topological_Data/TMData-%s_NS%d_NF%d.mat', Shuff_Method, n_sample, n_feature);
fprintf('Saving data in [%s]\n', sav_name);
save(sav_name, 'TM_Data', 'TM_Data_filtered', 'TM_Data_z', 'TM_Label', 'Pair_Info', 'TM_Name');
