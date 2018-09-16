clc;
clear;

%% Initialization
synet_fname = '../20_Prepare_SyNet_for_Cytoscape/Cyto_Input/SyNet_Top03544_Edge.tsv';
MAX_N_PAIR = inf;
rival_fname_lst = {
    'Basal_A', '../../Networks/StSpNet/A_core.txt'; ...
    'Basal_B', '../../Networks/StSpNet/‌‌‌‌‌‌B_core.txt'; ...
    'Luminal', '../../Networks/StSpNet/L_core.txt'; ...
    'Combined', '';
    };
n_rival = numel(rival_fname_lst(:,1));

%% Load SyNet
fid = fopen(synet_fname, 'r');
synet_pairs = textscan(fid, '%s%s%*s%*s', MAX_N_PAIR, 'Delimiter', '\t', 'ReturnOnError', 0, 'HeaderLine', 1);
if ~feof(fid), error('File is corrupted...'); end
n_pair = numel(synet_pairs{1});
SyNet_pair_lst = cell(n_pair, 1);
for pi=1:n_pair
    pair_item = sort({synet_pairs{1}{pi}, synet_pairs{2}{pi}});
    SyNet_pair_lst{pi} = sprintf('%s;%s', pair_item{1}, pair_item{2});
end
SyNet_pair_lst = unique(SyNet_pair_lst);
SyNet_gene_lst = unique([synet_pairs{:}]);

%% Loop over networks
ovl_info = table();
all_pairs_lst = {{}, {}};
for ni=1:n_rival
    fprintf('Loading %s network...\n', rival_fname_lst{ni, 1});
    if ~ strcmp(rival_fname_lst{ni, 1}, 'Combined')
        fid = fopen(rival_fname_lst{ni, 2}, 'r');
        riv_pairs = textscan(fid, '%s%*s%s', MAX_N_PAIR, 'Delimiter', ' ', 'ReturnOnError', 0);
        if ~feof(fid), error('File is corrupted...'); end
        all_pairs_lst{1} = [all_pairs_lst{1}; riv_pairs{1}];
        all_pairs_lst{2} = [all_pairs_lst{2}; riv_pairs{2}];
    else
        riv_pairs = all_pairs_lst;
    end
    n_pair = numel(riv_pairs{1});
    
    riv_pair_lst = cell(n_pair, 1);
    for pi=1:numel(riv_pairs{1})
        pair_item = sort({riv_pairs{1}{pi}, riv_pairs{2}{pi}});
        riv_pair_lst{pi} = sprintf('%s;%s', pair_item{1}, pair_item{2});
    end
    riv_pair_lst = unique(riv_pair_lst);
    riv_gene_lst = unique([riv_pairs{:}]);
    
    n_ovl_gene = numel(intersect(SyNet_gene_lst, riv_gene_lst));
    n_ovl_pair = numel(intersect(SyNet_pair_lst, riv_pair_lst));
    jak_ovl_gene = n_ovl_gene / numel(union(SyNet_gene_lst, riv_gene_lst));
    jak_ovl_pair = n_ovl_pair / numel(union(SyNet_pair_lst, riv_pair_lst));
    
    ovl_info = [ovl_info; ...
        cell2table({rival_fname_lst{ni, 1}, numel(riv_gene_lst), n_ovl_gene, numel(riv_pair_lst), n_ovl_pair, jak_ovl_gene, jak_ovl_pair}, ...
        'VariableNames', {'NetName', 'nGene', 'nGeneOvl', 'nPair', 'nPairOvl', 'Jak_Gene', 'Jak_Pair'})];
end
disp(ovl_info);

%% Save overlapping pairs
ovl_pairs_lst = intersect(SyNet_pair_lst, riv_pair_lst);
fid = fopen('./Output/Overlap_SyNet_SSBC.tsv', 'w');
for pi=1:numel(ovl_pairs_lst)
    pair_name = regexp(ovl_pairs_lst{pi}, ';', 'split');
    fprintf(fid, '%s\t%s\n', pair_name{1}, pair_name{2});
end
fclose(fid);