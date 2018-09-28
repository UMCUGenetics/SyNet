clc;
clear;

%% Initialization
code_path = '../11_Perform_LassoTypes/';
net_lst = {
    'SyNet-AvgSynACr-P50000'
    'STRING-P50000'
    };
n_net = numel(net_lst);
n_file = 140;

%% Loop over files
n_gene = zeros(n_file, n_net);
n_link = zeros(n_file, n_net);
for ni=1:n_net
    fn_lst = dir([code_path sprintf('Results_Files/DID_SyNet-SyNet_CVT01_*_%s_*_MSN-500_MTN-NetGL.mat', net_lst{ni})]);
    assert(n_file == numel(fn_lst));
    for fi=1:n_file
        fname = [code_path 'Results_Files/' fn_lst(fi).name];
        fprintf('Loading %s\n', fname);
        data = load(fname);
        n_gene(fi, ni) = numel(data.Gene_Name);
        
        SubNet_List = data.SubNet_List;
        n_sn = numel(SubNet_List);
        adj_net = zeros(n_gene(fi));
        for si=1:n_sn
            gene_set = SubNet_List{si};
            for gi=1:numel(gene_set)
                for gj=gi+1:numel(gene_set)
                    adj_net(gene_set(gi), gene_set(gj)) = 1;
                    adj_net(gene_set(gj), gene_set(gi)) = 1;
                end
            end
        end
        n_link(fi, ni) = sum(sum(triu(adj_net, 1)));
    end
end
