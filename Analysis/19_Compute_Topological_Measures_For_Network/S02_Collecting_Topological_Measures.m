clc;
clear;

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
n_pair = 10000;
n_lnk = 25000;

%% Load labels
SyNet_info = load('../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet_AvgSynACr.mat', 'PP_Info', 'NP_Info', 'Gene_Name');
Pair_Info = [
    SyNet_info.PP_Info(1:n_pair,:)  ones(n_pair, 1);
    SyNet_info.NP_Info(1:n_pair,:) -ones(n_pair, 1);
    ];
TM_Label = double(Pair_Info(:,16));

%% Collect data
net_lst = {'HumanInt' 'BioPlex','BioGRID','IntAct','STRING','HBBrain','HBKidney','HBOvary','HBLympNode','HBGland'}; % No 'AbsCorr'! as its used in SyNet
n_net = numel(net_lst);
tm_lst = {'ShortestPath' 'PageRank-FB0.95' 'PageRank-FB0.85' 'PageRank-FB0.75' 'PageRank-FB0.65' 'Eigenvector' 'Degree' 'Closeness' 'Betweenness'};
n_tm = numel(tm_lst);
TM_Data = zeros(n_pair*2, n_net*n_tm*2);
TM_Name = cell(n_net*n_tm*2, 1);
step = 1;
for ni=1:n_net
    for ti=1:n_tm
        data_name = sprintf('./Topological_Results/TM_%s_NP%06d_%s.mat', net_lst{ni}, n_lnk, tm_lst{ti});
        fprintf('Loading feature from [%s].\n', data_name);
        Data_Info = load(data_name);
        if ~isequal(SyNet_info.Gene_Name, Data_Info.Ref_GeneName), error(); end
        
        TM_Name{step}   = sprintf('%s-%s-A', net_lst{ni}, tm_lst{ti});
        TM_Data(:,step)   = Data_Info.Pair_AvgScore;
        TM_Name{step+1} = sprintf('%s-%s-D', net_lst{ni}, tm_lst{ti});
        TM_Data(:,step+1) = Data_Info.Pair_DifScore;
        step = step + 2;
    end
end
[n_sample, n_feature] = size(TM_Data);

%% Normalization
%imagesc(TM_Data);
TM_Data_filtered = TM_Data;
for fi=1:n_feature
    is_inf = TM_Data_filtered(:,fi)==inf;
    if sum(is_inf)>n_sample*0.5
        fprintf('[w] Warning: Large number of samples in [%s] are inf [%d/%d] (%0.2f%%).\n', TM_Name{fi}, sum(is_inf), n_sample, sum(is_inf)*100/n_sample);
    end
    TM_Data_filtered(is_inf,fi) = max(TM_Data_filtered(~is_inf,fi))*1.01;
end
if any(isnan(TM_Data_filtered(:))) || any(TM_Data_filtered(:)==-inf), error(); end
TM_Data_z = zscore(TM_Data_filtered);
%imagesc(TM_Data_z);

%% Saving Data
sav_name = sprintf('./Topological_Data/TMData_NS%d_NF%d.mat', n_sample, n_feature);
fprintf('Saving data in [%s]\n', sav_name);
save(sav_name, 'TM_Data', 'TM_Data_filtered', 'TM_Data_z', 'TM_Label', 'Pair_Info', 'TM_Name');
