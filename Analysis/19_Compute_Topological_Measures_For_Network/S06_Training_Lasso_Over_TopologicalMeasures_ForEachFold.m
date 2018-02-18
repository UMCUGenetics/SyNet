clc;
clear;
close all

%% Initialization
addpath('../_Utilities/');
addpath('../../../../Useful_Sample_Codes/ShowProgress');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
data_path = '../11_Perform_LassoTypes/Dataset_Files/';
IS_STRICT_CV = 1;
CLS_Name = 'Lasso';
n_Fold = 50;
% Shuff_Method = 'LnkShuff';
Shuff_Method = 'OneGRND';
n_pair = 7088;
n_study = 14;
n_rep = 10;

%% Load TM data
tmdata_ptr = dir('./Topological_Data/TMData-*');
TM_Data_z = [];
Pair_Info = [];
for ri=1:numel(tmdata_ptr)
    tmdata_name = ['./Topological_Data/' tmdata_ptr(ri).name];
    fprintf('Loading [%s] ...\n', tmdata_name);
    res_data = load(tmdata_name, 'TM_Data_z', 'TM_Label', 'TM_Name', 'Pair_Info');
    TM_Data_z = [TM_Data_z; res_data.TM_Data_z];
    Pair_Info = [Pair_Info; res_data.Pair_Info];
end

%% Main loop
for si=1:n_study
    for ri=1:n_rep
        %% Load dataset file
        data_ptr = sprintf('DID_CVT01_Si%02d-Ri%03d_AvgSynACr-P50000_*.mat', si, ri);
        data_info = dir([data_path data_ptr]);
        data_name = [data_path data_info.name];
        fprintf('Loading data from [%s]', data_name);
        data_info = load(data_name, 'DatasetTr');
        
        %% Get positive pairs
        Net_Adj = triu(data_info.DatasetTr.Net_Adj, 1);
        [Net_val, Net_SInd] = sort(Net_Adj(:), 'Descend');
        SyNet_PIndex = zeros(n_pair/2, 2);
        [SyNet_PIndex(:,1), SyNet_PIndex(:,2)] = ind2sub(size(Net_Adj), Net_SInd(1:n_pair/2));
        
        %% Collect the TM data
        has_TM = ismember(SyNet_PIndex, Pair_Info(:,1:2), 'rows');
        %% Evaluation of classifier
        Fold_Index = crossvalind('KFold', TM_Label, n_Fold);
        Fold_auc = zeros(n_Fold, 1);
        TM_PLabel = zeros(n_pair, 1);
        for fi=1:n_Fold
            iTr = Fold_Index~=fi;
            iTe = Fold_Index==fi;
            
            if IS_STRICT_CV
                Test_Index = Pair_Info(iTe,1:2);
                is_in = any(ismember(Pair_Info(:,1:2), Test_Index(:)), 2);
                iTr(is_in) = 0;
            end
            zTr = TM_Data_z(iTr, :);
            zTe = TM_Data_z(iTe, :);
            lTr = TM_Label(iTr);
            lTe = TM_Label(iTe);
            
            switch CLS_Name
                case 'Lasso'
                    lasso_opt = {'lassoType', 't', 'CV', 5, 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};
                    [B, fit] = lassoEx(zTr, lTr, lasso_opt{:});
                    opt_L = fit.IndexMinMSE;
                    pred_vec = zTe*B(:,opt_L);
                    [te_auc, optThresh] = getAUC(lTe, pred_vec, 50);
                    Fold_auc(fi) = te_auc * 100;
                    TM_PLabel(iTe) = (pred_vec >= optThresh)*2-1;
                otherwise
                    error();
            end
        end
        
        %% Save results
        sav_name = sprintf('./Saved_Pred/PredTM_CL-%s_SM-%s_NS-%d_NF-%d_ID-%06.0f.mat', CLS_Name, Shuff_Method, n_pair, n_Fold, rand*1e6);
        fprintf('Saving result of full TM pred in [%s]\n', sav_name);
        save(sav_name, 'TM_PLabel', 'TM_Label', 'Pair_Info');
    end
end