clc;
clear;

%% Initialization
met_lst = {'Avg', 'Std', 'DA2', 'Reg', 'PCA1', 'DPCA'};
net_lst = {'Random-NN05' 'Random-NN10' 'Random-NN20'};
res_path = './Result_Files/';
cv_ind = 50;
n_study = 14;
n_rep = 5;
n_epoch = 10000;
cmb_path = './Combined_AUC/';
[~, ~] = mkdir(cmb_path);

%% Main loop
for mi=1:numel(met_lst)
    for ni=1:numel(net_lst)
        res_ptr = sprintf([res_path 'RES_%s_%s_SyNet-SyNet_CVT%02d_*.mat'], met_lst{mi}, net_lst{ni}, cv_ind);
        res_info = dir(res_ptr);
        if numel(res_info)==0
            continue;
        end
        sav_name = sprintf([cmb_path 'CMB_%s_%s_CVT%02d.mat'], met_lst{mi}, net_lst{ni}, cv_ind);
        if exist(sav_name, 'file')
            fprintf('[i] Warning: Results are already collected for [%s], ignoring.\n', sav_name);
            continue;
        end
        
        fprintf('Getting results for [%s, %s] ... \n', met_lst{mi}, net_lst{ni});
        Te_AUC = zeros(n_epoch, n_rep, n_study, 'single');
        Ind_AUC = zeros(n_epoch, n_rep, n_study, 'single');
        for si=1:n_study
            for ri=1:n_rep
                res_name = sprintf([res_path 'RES_%s_%s_SyNet-SyNet_CVT%02d_Si%02d-Ri%03d.mat'], met_lst{mi}, net_lst{ni}, cv_ind, si, ri);
                %fprintf('\t\tLoading [%s]\n', res_name);
                res_info = load(res_name, 'Te_AUC', 'Ind_AUC');
                Te_AUC(:, ri, si) = res_info.Te_AUC;
                Ind_AUC(:, ri, si) = res_info.Ind_AUC;
            end
        end
        
        %% Saving results
        fprintf('Saving in [%s]\n', sav_name);
        save(sav_name, 'Te_AUC', 'Ind_AUC');
    end
end