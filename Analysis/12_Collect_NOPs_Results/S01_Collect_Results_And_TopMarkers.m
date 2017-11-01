% clc;
clear;

%% Initialization
result_path = '../11_Perform_LassoTypes/Results_Files/';
sav_path = './Collected_Results/';
% method_lst = {'TReg' 'TNMCAd' 'KNN0' 'KNN1' 'KNN3' 'KNN5' 'KNN7' 'TNMC'  'LExAG' 'TNMCAd' 'TLEx' ...
%     'DA2Lex' 'Lasso' 'GLasso' 'CFGLasso' 'GLasso2' 'GLasso7' 'GLasso10'
%     'GLasso20' 'TSVM-RBF' 'TKNN0' 'RnFrst' 'TRnFrst'
%     };
% net_lst = {
%     'None-G11748', 'Random-G00500', 'Random-P10000', ...
%     'AbsCorr-G00500', 'AbsCorr-P10000', 'ACr-P10000', ...
%     'AvgSyn-G00500', 'AvgSyn-P10000', 'AvgSynACr-G00500', 'AvgSynACr-P10000', ...
%     'HPRD-G00500', 'HPRD-P10000', 'I2D-G00500', 'I2D-P10000', 'KEGG-G00500', 'KEGG-P10000', 'MSigDB-G00500', 'MSigDB-P10000', ...
%     'STRING-G00500', 'STRING-P10000' ...
%     'HPRD-G11748' 'I2D-G11748' 'KEGG-G11748' 'STRING-G11748' 'MSigDB-G11748' ...
%     };
method_lst = {'KNN0' 'TKNN0' 'SVM-Lin' 'TSVM-Lin' 'SVM-RBF' 'TSVM-RBF' 'GLasso'}; % 
net_lst = {'None-G11748' 'ACr-G00500'};
feat_lst = [20 50 100 500 1000];
n_net = numel(net_lst);
n_met = numel(method_lst);
n_feat = numel(feat_lst);
cv_ind = 1;
n_study = 14;
n_rep = 10;
n_TopMarker = 500;
[~,~] = mkdir(sav_path);

%% Generate gene map
GE_Data = load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name', 'Study_Name');
n_gene = numel(GE_Data.Gene_Name);
global Gene_Map
Gene_Map = containers.Map(GE_Data.Gene_Name, 1:n_gene);

%% Main loop
for mi=1:n_met
    for ni=1:n_net
        for fi=1:n_feat
            sav_name = sprintf([sav_path 'MRK_CVT%02d_%s_%s_MSN-%03d.mat'], cv_ind, method_lst{mi}, net_lst{ni}, feat_lst(fi));
            if exist(sav_name, 'file')
                fprintf('[i] Results are already collected for [%s], ignoring ... \n', sav_name);
                continue;
            else
                res_ptr = sprintf('%sDID_CVT%02d_Si*-Ri*_%s_*_MSN-%03d_MTN-%s.mat', ...
                            result_path, cv_ind, net_lst{ni}, feat_lst(fi), method_lst{mi});
                file_info = dir(res_ptr);
                if numel(file_info)==0
                    continue;
                end
                fprintf('/// Collecting [%s]\n', sav_name);
                out_cmb = [];
                out_cmb.Marker_lst = [];
                out_cmb.Marker_SiRi = [];
                out_cmb.Marker_Score = [];
                out_cmb.AUC_mat = nan(n_study, n_rep);
                out_cmb.Gene_Map = Gene_Map;
                for si=1:n_study
                    for ri=1:n_rep
                        res_ptr = sprintf('%sDID_CVT%02d_Si%02d-Ri%03d_%s_*_MSN-%03d_MTN-%s.mat', ...
                            result_path, cv_ind, si, ri, net_lst{ni}, feat_lst(fi), method_lst{mi});
                        file_info = dir(res_ptr);
                        if numel(file_info)==0
                            fprintf('--- Missing [%60s]\n', res_ptr);
                            continue;
                        elseif numel(file_info)>1
                            fprintf('[i] Warning, multiple files found for [%60s]\n', res_ptr);
                            [~, sind] = sort([file_info(:).datenum], 'Descend');
                            file_info = file_info(sind);
                            for i=2:numel(file_info)
                                del_fname = [result_path file_info(i).name];
                                fprintf('Deleting [%s]\n', del_fname);
                                delete(del_fname);
                            end
                        end
                        fprintf('Reading from: [%s]\n', file_info(1).name);
                        res_data = load([result_path file_info(1).name]);
                        switch method_lst{mi}
                            case 'TNMC'
                                SubNet_Score = res_data.SubNet_Score;
                            case 'TNMCAd'
                                SubNet_Score = res_data.SubNet_Score(1:res_data.opt_K);
                                if ~isfield(out_cmb, 'Opt_K')
                                    out_cmb.Opt_K = res_data.opt_K;
                                else
                                    out_cmb.Opt_K(end+1,1) = res_data.opt_K;
                                end
                            case {'TLEx' 'TReg' 'DA2Lex'}
                                SubNet_Score = abs(res_data.B(:,res_data.fit.IndexMinMSE));
                            case 'Lasso'
                                SubNet_Score = abs(res_data.B(:,res_data.fit.IndexMinMSE));
                                res_data.SubNet_List = num2cell(1:numel(res_data.Gene_Name))';
                            case {'GLasso' 'GLasso2' 'GLasso7' 'GLasso10' 'GLasso20' 'CFGLasso' 'FERALAvgStdInt' 'FERALAvgStd' 'FERALInt'}
                                Group_Index = res_data.fit.Options.ind; %res_data.SubNet_GInd;
                                n_snet = size(Group_Index,2);
                                SubNet_Score = zeros(n_snet, 1);
                                for gi=1:n_snet
                                    set_ind = Group_Index(1,gi):Group_Index(2,gi);
                                    SubNet_Score(gi) = max(abs(res_data.B(set_ind, res_data.fit.IndexMinMSE)));
                                end
                            case 'LExAG'
                                SubNet_Score = abs(res_data.B(:,res_data.fit.IndexMinMSE));
                                res_data.SubNet_List = num2cell(1:n_gene)';
                                res_data.Gene_Name = GE_Data.Gene_Name;
                            case {'KNN0','KNN1','KNN3','KNN5','KNN7','TKNN0'}
                                n_mrk = numel(res_data.Gene_Name);
                                SubNet_Score = res_data.SubNet_Score(1:n_mrk);
                                res_data.SubNet_List = res_data.SubNet_List(1:n_mrk);
                                if ~isfield(out_cmb, 'Opt_K')
                                    out_cmb.Opt_K = res_data.opt_K;
                                else
                                    out_cmb.Opt_K(end+1,1) = res_data.opt_K;
                                end
                            case {'SVM-Lin' 'TSVM-Lin' 'SVM-RBF' 'TSVM-RBF'}
                                n_mrk = numel(res_data.Gene_Name);
                                SubNet_Score = res_data.SubNet_Score(1:n_mrk);
                                res_data.SubNet_List = res_data.SubNet_List(1:n_mrk);
                            case {'RnFrst' 'TRnFrst'}
                                n_mrk = numel(res_data.Gene_Name);
                                SubNet_Score = res_data.SubNet_Score(1:n_mrk);
                                res_data.SubNet_List = res_data.SubNet_List(1:n_mrk);
                            otherwise
                                error('Undefined method');
                        end
                        [SubNet_Score, SN_sid] = sort(SubNet_Score, 'Descend');
                        SN_sid(SubNet_Score==0) = [];
                        if numel(SN_sid)>n_TopMarker
                            SN_sid = SN_sid(1:n_TopMarker);
                        end
                        n_SN = min([numel(res_data.SubNet_List) numel(SN_sid)]);
                        Gene_Score = rand(n_gene, 1)*1e-10; % To make sure stability is not due to selection of all gene all the time
                        for gi=1:n_SN
                            g_ind = getGeneIndex(res_data.Gene_Name(res_data.SubNet_List{SN_sid(gi)}));
                            Gene_Score(g_ind) = Gene_Score(g_ind) + n_TopMarker - gi + 1;
                        end
                        [TopGene_Val, TopGene_Ind] = sort(Gene_Score, 'Descend');
                        TopGene_Ind(TopGene_Val<1) = 0;
                        out_cmb.Marker_lst(end+1, :) = TopGene_Ind(1:n_TopMarker)';
                        out_cmb.Marker_Score(end+1,:) = TopGene_Val(1:n_TopMarker)';
                        out_cmb.Marker_SiRi(end+1,:) = [si ri];
                        out_cmb.AUC_mat(si, ri) = res_data.te_auc;
                    end
                end
                if any(isnan(out_cmb.AUC_mat(:)))
                    fprintf('********* Data can not be saved .. NAN are found ... \n');
                else
                    save(sav_name, '-struct', 'out_cmb');
                end
            end
        end
    end
end

%% //////////////////// Functions
function g_ind = getGeneIndex(Gene_List)
global Gene_Map
g_ind = zeros(size(Gene_List));
for gi=1:numel(g_ind)
    g_ind(gi) = Gene_Map(Gene_List{gi});
end
end
