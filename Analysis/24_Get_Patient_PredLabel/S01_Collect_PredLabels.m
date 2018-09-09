% clc;
clear;

%% Initialization
result_path = '../11_Perform_LassoTypes/Results_Files/';
sav_fname = './Predicted_PatientLabels.xlsx';
method_lst = {'LExAG', 'NetGL'}; % , 'NetLasso'
net_lst = {'AvgSynACr-P50000'};
n_net = numel(net_lst);
n_met = numel(method_lst);
cv_ind = 1;
n_study = 14;
n_rep = 10;
auc_mat = nan(n_met, n_rep, n_study);
GEPath = '../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat';

%% Load Expression data
fprintf('Loading train expression data from [%s] ...\n', GEPath);
expr_info = load(GEPath, 'Patient_Label', 'Patient_Info');
assert(isequal(expr_info.Patient_Info(:, 'Prognostic_Status').Variables, expr_info.Patient_Label));
n_pat = length(expr_info.Patient_Label);
Pred_tbl = expr_info.Patient_Info(:, {'PatientID', 'Subtype', 'SurvivalTime', 'Prognostic_Status'});

%% Main loop
for mi=1:n_met
    for ni=1:n_net
        for ri=1:n_rep
%             col_sid = strrep(sprintf('%s_%s_Rep%02d', method_lst{mi}, net_lst{ni}, ri), '-', '_');
            col_sid = strrep(sprintf('%s_Rep%02d', method_lst{mi}, ri), '-', '_');
            Pred_tbl = addvars(Pred_tbl, nan(n_pat, 1), 'NewVariableNames', col_sid);
            for si=1:n_study
                res_ptr = sprintf('%sDID_CVT%02d_Si%02d-Ri%03d_%s_*_MSN-500_MTN-%s.mat', ...
                    result_path, cv_ind, si, ri, net_lst{ni}, method_lst{mi});
                file_info = dir(res_ptr);
                if numel(file_info)~=1
                    error('--- Missing [%60s]\n', res_ptr);
                end
                fprintf('Reading from: [%s]\n', file_info(1).name);
                res_data = load([result_path file_info(1).name]);
                
                %% Process results
                assert(all(isnan(Pred_tbl(res_data.UsedTeSamples, col_sid).Variables)));
                [te_auc, optThresh] = getAUC(res_data.TeLbl, res_data.predTeLbl, 50);
                assert(abs(te_auc - res_data.te_auc) < 1e-10);
                Pred_tbl(res_data.UsedTeSamples, col_sid) = array2table(res_data.predTeLbl > optThresh);
                
                auc_mat(mi, ri, si) = te_auc;
            end
        end
    end
end

%% Export collected information
writetable(Pred_tbl, sav_fname);


