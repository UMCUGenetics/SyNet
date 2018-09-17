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
out_name = './Output/Overlap_SyNet_SSBC';

%% Load SyNet
fid = fopen(synet_fname, 'r');
file_cell = textscan(fid, '%s%s%*s%*s', MAX_N_PAIR, 'Delimiter', '\t', 'ReturnOnError', 0, 'HeaderLine', 1);
if ~feof(fid), error('File is corrupted...'); end
n_pair = numel(file_cell{1});
synet_sid_plst = cell(n_pair, 1);
for pi=1:n_pair
    pair_item = sort({file_cell{1}{pi}, file_cell{2}{pi}});
    synet_sid_plst{pi} = sprintf('%s;%s', pair_item{1}, pair_item{2});
end
synet_sid_plst = unique(synet_sid_plst);
synet_gene_lst = unique([file_cell{:}]);

%% Loop over networks
ovl_info = table();
all_pairs_lst = {{}, {}};
for ni=1:n_rival
    fprintf('Loading %s network...\n', rival_fname_lst{ni, 1});
    if ~ strcmp(rival_fname_lst{ni, 1}, 'Combined')
        fid = fopen(rival_fname_lst{ni, 2}, 'r');
        file_cell = textscan(fid, '%s%*s%s', MAX_N_PAIR, 'Delimiter', ' ', 'ReturnOnError', 0);
        if ~feof(fid), error('File is corrupted...'); end
        all_pairs_lst{1} = [all_pairs_lst{1}; file_cell{1}];
        all_pairs_lst{2} = [all_pairs_lst{2}; file_cell{2}];
    else
        file_cell = all_pairs_lst;
    end
    n_pair = numel(file_cell{1});
    
    riv_sid_plst = cell(n_pair, 1);
    for pi=1:numel(file_cell{1})
        pair_item = sort({file_cell{1}{pi}, file_cell{2}{pi}});
        riv_sid_plst{pi} = sprintf('%s;%s', pair_item{1}, pair_item{2});
    end
    riv_sid_plst = unique(riv_sid_plst);
    riv_gene_lst = unique([file_cell{:}]);
    
    n_ovl_gene = numel(intersect(synet_gene_lst, riv_gene_lst));
    n_ovl_pair = numel(intersect(synet_sid_plst, riv_sid_plst));
    jak_ovl_gene = n_ovl_gene / numel(union(synet_gene_lst, riv_gene_lst));
    jak_ovl_pair = n_ovl_pair / numel(union(synet_sid_plst, riv_sid_plst));
    
    ovl_info = [ovl_info; ...
        cell2table({rival_fname_lst{ni, 1}, numel(riv_gene_lst), n_ovl_gene, jak_ovl_gene, numel(riv_sid_plst), n_ovl_pair, jak_ovl_pair}, ...
        'VariableNames', {'NetworkName', 'n_gene', 'n_gene_overlap', 'gene_ovl_Jaccard', 'n_pair', 'n_pair_overlap', 'pair_ovl_Jaccard'})];
end
disp(ovl_info);
writetable(ovl_info, [out_name '.xlsx']);

%% Save overlapping pairs
ovl_pairs_lst = intersect(synet_sid_plst, riv_sid_plst);
fid = fopen([out_name '.tsv'], 'w');
fprintf(fid, 'Source\tTarget\tType\tWeight\n');
for pi=1:numel(ovl_pairs_lst)
    pair_name = regexp(ovl_pairs_lst{pi}, ';', 'split');
    fprintf(fid, '%s\t%s\tUndirected\t1\n', pair_name{1}, pair_name{2});
end
fclose(fid);

