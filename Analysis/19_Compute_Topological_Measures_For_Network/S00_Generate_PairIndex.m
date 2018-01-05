clc;
clear;

%% Inialization
Ref_Name = 'AvgSynACr';
MAX_SyNet_Pair = 20000;
Shuff_Method = 'LnkShuff';

%% Load SyNet pairs
SyNet_name = sprintf('../01_Pairwise_Evaluation_of_Genes/Top_Pairs/TopP_SyNet_%s.mat', Ref_Name);
fprintf('Loading [%s] as reference network ...\n', Ref_Name);
SyNet_info = load(SyNet_name, 'PP_Info', 'NP_Info', 'Gene_Name');
n_gene = numel(SyNet_info.Gene_Name);
PP_Info = SyNet_info.PP_Info(1:MAX_SyNet_Pair, :);

%% Generate Pair index
switch Shuff_Method
    case 'LnkShuff'
        row_rind = randperm(MAX_SyNet_Pair);
        NP_Info = [PP_Info(:,1) PP_Info(row_rind,2)];
        if ~isequal(sort(NP_Info), sort(PP_Info(:,1:2))), error(); end
    case 'OneGRND'
        NP_Info = zeros(MAX_SyNet_Pair, 2);
        for pi=1:MAX_SyNet_Pair
            col_rind = randperm(2,1);
            NP_Info(pi, :) = [PP_Info(pi, col_rind) randperm(n_gene,1)];
        end
end

%% Saving indices
PI_Name = sprintf('./Topological_Data/PairInfo-%s_%s_MP%06d.mat', Shuff_Method, Ref_Name, MAX_SyNet_Pair);
save(PI_Name, 'SyNet_name', 'SyNet_info', 'PP_Info', 'NP_Info');
