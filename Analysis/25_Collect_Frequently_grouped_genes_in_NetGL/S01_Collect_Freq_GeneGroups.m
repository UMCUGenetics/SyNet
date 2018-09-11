% clc;
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
n_top = 500;
out_fname = sprinf([sav_path 'Grp_%s_%s.mat'], net_name, method_name);

%% Generate gene map
GE_Data = load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name', 'Study_Name');
n_gene = numel(GE_Data.Gene_Name);
Gene_Map = containers.Map(GE_Data.Gene_Name, 1:n_gene);

%% Main loop
Nei_mat = zeros(n_gene, uint8);
for ri=1:n_rep
    for fi=1:n_study
        sav_name = sprintf([sav_path 'Grp_CVT%02d_%s_%s_MSN-%03d.mat'], cv_ind, method_name, net_name, MAX_N_SUBNET);
        if exist(sav_name, 'file')
            fprintf('[i] Results are already collected for [%s], ignoring ... \n', sav_name);
            continue;
        else
            res_ptr = sprintf('%sDID_CVT%02d_Si%02d-Ri%03d_%s_*_MSN-%03d_MTN-%s.mat', ...
                result_path, cv_ind, fi, ri, net_name, MAX_N_SUBNET, method_name);
            file_info = dir(res_ptr);
            if numel(file_info)==0
                fprintf('[w] No results found for [%s] ...\n', res_ptr);
                continue;
            end
            res_fname = [result_path file_info(1).name];
            fprintf('/// Loading [%s]\n', res_fname);
            res_data = load(res_fname);
            n_SubNet = length(res_data.SubNet_List);
            for si=1:n_SubNet
                gene_set = res_data.SubNet_List{si};
                grp_size = length(gene_set);
                gene_index = zeros(grp_size, 1);
                for gi=1:grp_size
                    gname = res_data.Gene_Name{gene_set(gi)};
                    gene_index(gi) = Gene_Map(gname);
                end
                for gi=1:grp_size
                    for gj=1:grp_size
                        Nei_mat(gene_index(gi), gene_index(gj)) = Nei_mat(gene_index(gi), gene_index(gj)) + 1;
                        Nei_mat(gene_index(gj), gene_index(gi)) = Nei_mat(gene_index(gj), gene_index(gi)) + 1;
                    end
                end
            end
        end
    end
end

%% Identify the freq neighbors
assert(issymmetric(Nei_mat));
Nei_mat = triu(Nei_mat, 1);
sorted_index = sort(Nei_mat(:), 'descend');
nei_freq = zeros(n_top, 1);
nei_plst = zeros(n_top, 2);
for ti=1:n_top
    p_ind = ind2sub(sorted_index, [n_gene, n_gene]);
    nei_freq(ti) = Nei_mat(p_ind(1), p_ind(2));
    nei_plst(ti, :) = p_ind;
end

%% Saving the collected results
save(out_fname, 'Gene_Map', 'nei_freq', 'nei_plst');

