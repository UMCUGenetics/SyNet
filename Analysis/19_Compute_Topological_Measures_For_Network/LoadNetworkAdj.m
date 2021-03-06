function net_info = LoadNetworkAdj(net_name, net_opt)

%% Load network data
net_info.net_name = net_name;
net_info.MAX_N_PAIR = net_opt.MAX_N_PAIR;
switch 1
    case ismember(net_name, {'AbsCorr'})
        ge_data = load(net_opt.GE_Path, 'Gene_Expression', 'Gene_Name');
        Net_Adj = abs(corr(zscore(ge_data.Gene_Expression), 'Type', 'Spearman'));
        Net_Adj(Net_Adj<0.5) = 0;
        Gene_Name = ge_data.Gene_Name;
        net_info.net_path = net_opt.GE_Path;
        clear ge_data
    case ismember(net_name, {'KEGG', 'MSigDB'})
        net_info.net_path = getPath(net_name);
        GSet_lst = regexp(fileread(net_info.net_path), '\n', 'split')';
        if strcmp(GSet_lst{end},''), GSet_lst(end)=[]; end
        n_gset = numel(GSet_lst);
        fprintf('Loading [%d] gene sets from [%s] ...\n', n_gset, net_info.net_path);
        for si=1:n_gset
            GSet_lst{si} = regexp(GSet_lst{si}, '\t', 'split');
        end
        Gene_Name = unique([GSet_lst{:}])';
        n_gene = numel(Gene_Name);
        fprintf('[i] Network contains [%d] genes before filtering.\n', n_gene);
        
        GMap = containers.Map(Gene_Name, 1:n_gene);
        Net_Adj = zeros(n_gene, 'single');
        fprintf('Adding genes to network: ');
        for si=1:n_gset
            showprogress(si, n_gset, 20);
            grp_size = numel(GSet_lst{si});
            g_ind = zeros(grp_size, 1);
            for gi=1:grp_size
                if GMap.isKey(GSet_lst{si}{gi})
                    g_ind(gi) = GMap(GSet_lst{si}{gi});
                end
            end
            g_ind = nonzeros(g_ind);
            Net_Adj(g_ind, g_ind) = rand(numel(g_ind));
        end
        Net_Adj = max(Net_Adj, Net_Adj');
        clear GSet_lst GMap
    case any(~cellfun('isempty', regexp(net_name, {'HB.*' 'STRING','STRINGnShuff','HPRD','I2D','IntAct','HumanInt','BioPlex','BioGRID'})))
        net_info.net_path = getPath(net_name);
        fid = fopen(net_info.net_path, 'r');
        Header_lst = regexp(fgetl(fid), '\t', 'split');
        if numel(Header_lst)==2
            fprintf('No link weight exists. Selection of [%d] links from begining of file.\n', net_info.MAX_N_PAIR*5);
            net_cell = textscan(fid,   '%s%s', net_info.MAX_N_PAIR*5, 'Delimiter', '\t', 'ReturnOnError', 0);
        else
            fprintf('Link weights are found. Selecting of [%d] links from top weighted interactions.\n', net_info.MAX_N_PAIR*5);
            net_cell = textscan(fid, '%s%s%f', net_info.MAX_N_PAIR*5, 'Delimiter', '\t', 'ReturnOnError', 0);
            if ~issorted(net_cell{3}, 'Descend'), error('Given network file is not sorted.'); end
        end
        fclose(fid);
        Gene_Name = unique(vertcat(net_cell{1:2}));
        n_gene = numel(Gene_Name);
        GMap = containers.Map(Gene_Name, 1:n_gene);
        Net_Adj = zeros(n_gene, 'single');
        n_lnk = numel(net_cell{1});
        fprintf('[%d] genes and [%d] links are loaded from [%s]\n', n_gene, n_lnk, net_name);
        fprintf('Forming the Adj matrix ...\n');
        for li=1:n_lnk
            showprogress(li, n_lnk);
            if GMap.isKey(net_cell{1}{li}) && GMap.isKey(net_cell{2}{li})
                gi = GMap(net_cell{1}{li});
                gj = GMap(net_cell{2}{li});
                Net_Adj(gi, gj) = n_lnk - li + 1;
                Net_Adj(gj, gi) = n_lnk - li + 1;
            end
        end
        clear net_cell GMap
    otherwise
        error('Unknown method.');
end
if ~issymmetric(Net_Adj), fprintf('Warning: Raw Adj matrix is not symetric.\n'); end
Net_Adj = max(Net_Adj, Net_Adj');
Net_Adj(1:size(Net_Adj,1)+1:end) = 0;

%% Unifying gene list
if isfield(net_opt, 'PreferredGenes')
    ReferenceGenes = intersect(Gene_Name, net_opt.PreferredGenes);
    Ind_List = UnifyGeneList(ReferenceGenes, Gene_Name);
    Net_Adj = Net_Adj(Ind_List, Ind_List);
    Gene_Name = Gene_Name(Ind_List);
    if ~isequal(ReferenceGenes, Gene_Name), error(); end
    fprintf('[%d] genes are left after filtering according to provided list.\n', numel(Gene_Name));
end

%% Link filtering
if isfield(net_info, 'MAX_N_PAIR')
    fprintf('[i] Selecting top %d interactions: ', net_info.MAX_N_PAIR);
    Tmp_Adj = triu(Net_Adj, 1);
    [~, scr_ind] = sort(Tmp_Adj(:), 'Descend');
    clear Tmp_Adj
    if numel(unique(Net_Adj(scr_ind(1:net_info.MAX_N_PAIR))))<=1
        error('Identical weights are found. This function is not implemented for that. Shuffle the links before removing');
    end
    scr_ind(1:net_info.MAX_N_PAIR) = [];
    net_info.Net_Threshold = Net_Adj(scr_ind(1));
    Net_Adj(scr_ind) = 0;
    Net_Adj = max(Net_Adj, Net_Adj');
    if net_info.Net_Threshold<=0
        fprintf('Warning: Identified threshhold is [%d]. Fewer pairs are selected than requested [%d]\n', net_info.Net_Threshold, net_info.MAX_N_PAIR);
    end
    clear scr_val scr_ind
    fprintf('[%d] genes and [%d] links are left in the network.\n', numel(Gene_Name), numel(nonzeros(triu(Net_Adj))));
end

%% Node filtering
fprintf('[i] Filtering nodes with no iteractions ...\n');
del_ind = sum(Net_Adj~=0,1)==0;
Net_Adj(del_ind, :) = [];
Net_Adj(:, del_ind) = [];
Gene_Name(del_ind) = [];

%% Storing
if ~issymmetric(Net_Adj), error('Adj Matrix is not symetric.\n'); end
net_info.Net_Adj = single(Net_Adj);
net_info.Gene_Name = Gene_Name;
net_info.N_PAIR = numel(nonzeros(triu(Net_Adj)));
fprintf('[%d] genes and [%d] links are left in the network.\n', numel(Gene_Name), net_info.N_PAIR);
end

