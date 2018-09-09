% clc;
clear;

%% Initialization
result_path = '../11_Perform_LassoTypes/Results_Files/';
sav_fname = './predicted_patient_labels.tsv';
method_lst = {'NetLasso', 'NetGL'}; %  'NetLasso' 'NetGL' 'CvGL' 'TMGL' 'Lasso' 'GLasso' 'CFGLasso' 'GLasso5'
net_lst = {'AvgSynACr-P50000'};
n_net = numel(net_lst);
n_met = numel(method_lst);
cv_ind = 1;
n_study = 14;
n_rep = 10;

%% Main loop
for mi=1:n_met
    for ni=1:n_net
        for si=1:n_study
            for ri=1:n_rep
                res_ptr = sprintf('%sDID_CVT%02d_Si%02d-Ri%03d_%s_*_MSN-500_MTN-%s.mat', ...
                    result_path, cv_ind, si, ri, net_lst{ni}, method_lst{mi});
                file_info = dir(res_ptr);
                if numel(file_info)~=1
                    error('--- Missing [%60s]\n', res_ptr);
                end
                fprintf('Reading from: [%s]\n', file_info(1).name);
                res_data = load([result_path file_info(1).name]);
                
                %% Process results
                res_data.fit.
            end
        end
    end
end



