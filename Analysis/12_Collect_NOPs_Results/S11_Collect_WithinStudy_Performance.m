clc;
clear;

%% Initialization
result_path = '../11_Perform_LassoTypes/Results_Files/';
sav_path = './Collected_Results/';
Net_lst = {'AvgSynACr-P01000','AvgSynACr-P05000','AvgSynACr-P10000','AvgSynACr-P20000','AvgSynACr-P50000' ...
    'AvgSynACr-P100000','AvgSynACr-P150000','AvgSynACr-P200000','AvgSynACr-P500000'...
    };
Met_lst = {'GLasso2' 'GLasso3' 'GLasso5' 'GLasso7' 'GLasso10' 'GLasso15' 'GLasso20'};
n_net = numel(Net_lst);
n_met = numel(Met_lst);
n_rep = 5;
study_index = 14;

%% Main loop
AUC_mat = nan(n_met, n_net, n_rep);
for ni=1:n_net
    for mi=1:n_met
        for ri=1:n_rep
            res_ptr = sprintf([result_path 'DID_CVT51_Si%02d-Ri%03d_%s_*_MSN-500_MTN-%s.mat'], study_index, ri, Net_lst{ni}, Met_lst{mi});
            res_info = dir(res_ptr);
            if numel(res_info)>1
                fprintf('[w] Multiple results found.\n%s\n', strjoin({res_info.name}, '\n'));
                del_name = [result_path res_info(1).name];
                fprintf('Deleting [%s] ...\n', del_name);
                delete(del_name);
                res_info(1) = [];
            end
            if numel(res_info)==1
                fprintf('Loading [%s]\n', res_info.name);
                res_info = load([result_path res_info.name]);
                AUC_mat(mi, ni, ri) = res_info.te_auc;
            else
                fprintf('[w] Warning: File [%s] not found ...\n', res_ptr);
            end
        end
    end
end

%% Save
out_name = sprintf([sav_path 'WithinStudy_Si%02d_Results.mat'], study_index);
fprintf('Saving results in [%s] ...\n', out_name);
save(out_name, 'AUC_mat', 'Net_lst', 'Met_lst');

return
%% Plotting results
clear
close all
res_info = load('./Collected_Results/WithinStudy_Si14_Results.mat');
imagesc(mean(res_info.AUC_mat, 3));
n_net = numel(res_info.Net_lst);
n_met = numel(res_info.Met_lst);
set(gca, 'XTick', 1:n_net, 'XTickLabel', res_info.Net_lst, 'XTickLabelRotation', 30, ...
'YTick', 1:n_met, 'YTickLabel', res_info.Met_lst, 'Clim', [0.68 0.80]);
colorbar();

