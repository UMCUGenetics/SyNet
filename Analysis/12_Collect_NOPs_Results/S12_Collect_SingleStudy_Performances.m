clc;
clear;

%% Initialization
result_path = '../11_Perform_LassoTypes/Results_Files/';
sav_path = './Collected_Results/';
Net_lst = {'AvgSynACr-P01000','AvgSynACr-P05000','AvgSynACr-P10000','AvgSynACr-P20000','AvgSynACr-P50000' ...
    'AvgSynACr-P100000','AvgSynACr-P150000','AvgSynACr-P200000','AvgSynACr-P500000','AvgSynACr-P1000000'...
    };
Met_lst = {'GLasso2' 'GLasso3' 'GLasso5' 'GLasso7' 'GLasso10' 'GLasso15' 'GLasso20'};
n_net = numel(Net_lst);
n_met = numel(Met_lst);
n_rep = 5;
study_index = 14;
CV_ID = 01;

%% Main loop
AUC_mat = nan(n_met, n_net, n_rep);
for ni=1:n_net
    for mi=1:n_met
        for ri=1:n_rep
            res_ptr = sprintf([result_path 'DID_CVT%02d_Si%02d-Ri%03d_%s_*_MSN-500_MTN-%s.mat'], CV_ID, study_index, ri, Net_lst{ni}, Met_lst{mi});
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
                fprintf('[w] Warning: Non-unique or missing file for [%s] pattern ...\n', res_ptr);
            end
        end
    end
end

%% Save
out_name = sprintf([sav_path 'SingleStudy_CID%02d_Si%02d_Results.mat'], CV_ID, study_index);
fprintf('Saving results in [%s] ...\n', out_name);
save(out_name, 'AUC_mat', 'Net_lst', 'Met_lst');

return

%% Plotting results
clear
close all
GE_Info = load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Study_Name');
figure();
study_index = 14;
CV_ID = 52;
res_name = sprintf('./Collected_Results/SingleStudy_CID%02d_Si%02d_Results.mat', CV_ID, study_index);
fprintf('Loading [%s]\n', res_name);
res_info = load(res_name);
imagesc(mean(res_info.AUC_mat, 3)*100);
n_net = numel(res_info.Net_lst);
n_met = numel(res_info.Met_lst);
set(gca, 'XTick', 1:n_net, 'XTickLabel', res_info.Net_lst, 'XTickLabelRotation', 20, ...
'YTick', 1:n_met, 'YTickLabel', res_info.Met_lst, 'FontWeight', 'Bold'); %
% xlim([0.5 4.5]);
clrmap_h = colormap(flipud(summer(10)));
title(sprintf('Performance of Group Lasso for [%s] study on CV %d', GE_Info.Study_Name{study_index}, CV_ID));
clrbar_h = colorbar();
ylabel(clrbar_h, 'Are under the curve', 'FontWeight', 'Bold');
xlabel('# Top pairs in SyNet', 'FontWeight', 'Bold');
ylabel('Gene set size', 'FontWeight', 'Bold');

%% Saving
output_name = sprintf('./Plots/S12_GLasso_SingleStudyPerf_CID%02d_Si%02d.pdf', CV_ID, study_index);
set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [20 15], 'PaperPosition', [0 0 20 15]);
print('-dpdf', '-r300', output_name);
