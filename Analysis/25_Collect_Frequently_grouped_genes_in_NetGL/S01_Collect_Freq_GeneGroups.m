clc;
clear;

%% Initialization
result_path = '../11_Perform_LassoTypes/Results_Files_afterAddingPredLabel_toResFiles/';
sav_path = './Collected_Results/';
method_name = 'NetGL';
net_name = 'SyNet-AvgSynACr-P50000';
MAX_N_SUBNET = 500;
cv_ind = 1;
n_study = 14;
n_rep = 10;
n_top = 10000;
out_fname = sprintf([sav_path 'Grp_CVT%02d_%s_%s_MSN-%03d_NT%d'], cv_ind, net_name, method_name, MAX_N_SUBNET, n_top);

%% Generate gene map
GE_Data = load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name', 'Study_Name');
n_gene = numel(GE_Data.Gene_Name);
assert(n_gene == numel(unique(GE_Data.Gene_Name)));
Gene_Map = containers.Map(GE_Data.Gene_Name, 1:n_gene);

%% Main loop
nei_mat = zeros(n_gene, 'uint16');
for ri=1:n_rep
    for si=1:n_study
        res_ptr = sprintf([result_path 'DID_CVT%02d_Si%02d-Ri%03d_%s_*_MSN-%03d_MTN-%s.mat'], ...
            cv_ind, si, ri, net_name, MAX_N_SUBNET, method_name);
        file_info = dir(res_ptr);
        if numel(file_info)~=1
            error('[w] No results found for [%s] ...\n', res_ptr);
        end
        res_fname = [result_path file_info(1).name];
        fprintf('/// Loading [%s]\n', res_fname);
        res_data = load(res_fname);
        n_SubNet = length(res_data.SubNet_List);
        for ni=1:n_SubNet
            run_gindex_lst = res_data.SubNet_List{ni};
            grp_size = length(run_gindex_lst);
            global_gindex_lst = zeros(grp_size, 1);
            for gi=1:grp_size
                gname = res_data.Gene_Name{run_gindex_lst(gi)};
                global_gindex_lst(gi) = Gene_Map(gname);
            end
            for gi=1:grp_size
                for gj=1:grp_size
                    nei_mat(global_gindex_lst(gi), global_gindex_lst(gj)) = nei_mat(global_gindex_lst(gi), global_gindex_lst(gj)) + 1;
                    nei_mat(global_gindex_lst(gj), global_gindex_lst(gi)) = nei_mat(global_gindex_lst(gj), global_gindex_lst(gi)) + 1;
                end
            end
        end
    end
end

%% Identify the freq neighbors
assert(issymmetric(double(nei_mat)));
nei_mat = triu(nei_mat, 1);
[sorted_freq, sorted_index] = sort(nei_mat(:), 'descend');
nei_freq = zeros(n_top, 1);
nei_plst = zeros(n_top, 2);
nei_gname = cell(n_top, 2);
for ti=1:n_top
    [p_i, p_j] = ind2sub([n_gene, n_gene], sorted_index(ti));
    assert(sorted_freq(ti) == nei_mat(p_i, p_j));
    nei_freq(ti) = nei_mat(p_i, p_j);
    nei_plst(ti, :) = [p_i, p_j];
    nei_gname(ti, :) = {GE_Data.Gene_Name{p_i}, GE_Data.Gene_Name{p_j}};
end
fprintf('Total of %d pairs and %d unique genes are selected.\n', n_top, length(unique(nei_plst)));

%% Saving the collected results
save([out_fname '.mat'], 'Gene_Map', 'nei_freq', 'nei_plst', 'nei_gname');

fid = fopen([out_fname '.tsv'], 'w');
fprintf(fid, 'Source\tTarget\tWeight\tType\n');
for ti=1:n_top
    fprintf(fid, '%s\t%s\t%d\tUndirected\n', ...
    nei_gname{ti, 1}, nei_gname{ti, 2}, nei_freq(ti));
end
fclose(fid);

