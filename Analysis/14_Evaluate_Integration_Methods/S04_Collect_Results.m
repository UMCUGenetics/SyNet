clc;
clear;

%% Initialization
res_path = './Result_Files/';
cmb_path = './Combined_AUC/';
[~, ~] = mkdir(cmb_path);
met_lst = {'Reg'}; % 'Avg', 'Std', 'DA2', 'Reg', 'PCA1', 'DPCA', 'DA2NoRem', 'Rnd'
net_lst = {'STRING' 'STRING-Shf'};
gset_lst = [1 2 5 7 10 20 50];
study_lst = [14];
n_study = numel(study_lst);
cv_ind = 51;
n_rep = 5;

%% Main loop
for gi=1:numel(gset_lst)
    for mi=1:numel(met_lst)
        for ni=1:numel(net_lst)
            res_ptr = sprintf([res_path 'RES_%s_%s_NN%02d_SyNet-SyNet_CVT%02d_*.mat'], met_lst{mi}, net_lst{ni}, gset_lst(gi), cv_ind);
            res_info = dir(res_ptr);
            if numel(res_info)==0
                continue;
            end
            sav_name = sprintf([cmb_path 'CMB_%s_%s_NN%02d_CVT%02d.mat'], met_lst{mi}, net_lst{ni}, gset_lst(gi), cv_ind);
            if exist(sav_name, 'file')
                fprintf('[i] Warning: Results are already collected for [%s], ignoring.\n', sav_name);
                continue;
            end
            
            fprintf('Getting results for [%s, %s, %02d] ... \n', met_lst{mi}, net_lst{ni}, gset_lst(gi));
            for si=1:n_study
                for ri=1:n_rep
                    res_name = sprintf([res_path 'RES_%s_%s_NN%02d_SyNet-SyNet_CVT%02d_Si%02d-Ri%03d.mat'], met_lst{mi}, net_lst{ni}, gset_lst(gi), cv_ind, study_lst(si), ri);
                    %fprintf('\t\tLoading [%s]\n', res_name);
                    res_info = load(res_name, 'Cmb_TeAUC', 'Ind_TeAUC');
                    if si==1 && ri==1
                        n_epoch = numel(res_info.Cmb_TeAUC);
                        Cmb_TeAUC = zeros(n_epoch, n_rep, n_study, 'single');
                        Ind_TeAUC = zeros(n_epoch, n_rep, n_study, 'single');
                    end
                    Cmb_TeAUC(:, ri, si) = res_info.Cmb_TeAUC;
                    Ind_TeAUC(:, ri, si) = res_info.Ind_TeAUC;
                end
            end
            
            %% Saving results
            fprintf('Saving in [%s]\n', sav_name);
            save(sav_name, 'Cmb_TeAUC', 'Ind_TeAUC');
        end
    end
end