function S00_Main_Code(Target_Study, Target_Repeat, method_lst, net_lst, MAX_N_SUBNET)
%% Run
%{
for ri in `seq 1 10`; do
for si in `seq 1 14`; do
PARAM="$si,$ri,{'NetLasso','NetGL'},{'AvgSynACr-P10000'}"; sbatch --exclude=maxwell --job-name=NE-$PARAM --output=Logs/NE-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=10GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S00_Main_Code "$PARAM";
done;
read -p "`date`: $PARAM. Press a key" -t 1800
done

UMC: PARAM="$si,$ri,{'TAgNMC','TNMC','TLEx','TAgLEx'},{'Random-T00010'},10"; qsub -N "NE-$PARAM" -l h_rt=24:00:00 -l h_vmem=5G ~/bulk/env/run_Matlab.sh S00_Main_Code "$PARAM";
%}

%% ####
if ismac || ispc
    fprintf('*** Warning!: Running on debug mode.\n');
    Target_Study = 14;
    Target_Repeat = 2;
    method_lst = {'NetLasso', 'NetGL'};
    net_lst = {'ACrNShuff-P25000'};
    MAX_N_SUBNET = 500;
end

%% Initialization
clc
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
addpath('../_Utilities/');
dataset_path = './Dataset_Files/';
result_path = './Results_Files/';
opt_info.lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};
if ~exist('MAX_N_SUBNET', 'var')
    opt_info.MAX_N_SUBNET = 500;
else
    opt_info.MAX_N_SUBNET = MAX_N_SUBNET;
end

n_net = numel(net_lst);
n_meth = numel(method_lst);

%% Main Loop
fprintf([repmat('/',1,20) ' Start of main loop ' repmat('/',1,20) '\n']);
cv_id = sprintf('CVT01_Si%02d-Ri%03d', Target_Study, Target_Repeat);
fprintf('[i] CV ID is: %s\n', cv_id);
fprintf('[i] Method list is: %s\n', strjoin(method_lst, ', '));
fprintf('[i] Network list is: %s\n', strjoin(net_lst, ', '));
fprintf([repmat('/',1,60) '\n']);

te_auc = nan(n_net, n_meth);
for ni=1:n_net
    fprintf(['****** Network Gen [%s] ' repmat('*',1,40) '\n'], net_lst{ni});
    S02_GenerateDataset(cv_id, net_lst{ni});
    
    %% Load Dataset info
    ds_id = sprintf('%s_%s', cv_id, net_lst{ni});
    dataset_list = dir([dataset_path 'DID_' ds_id '_*.mat']);
    if numel(dataset_list)~=1, error('Missing or duplicate dataset found.. \n[%s]\n', strjoin({dataset_list.name}, ', ')); end
    dataset_name = [dataset_path dataset_list(1).name];
    fprintf('Loading dataset [%s] ...\n', dataset_name);
    dataset_info = load(dataset_name);
    dataset_info.dataset_name = dataset_name;
    fprintf('Dataset has Train: [%d x %d], Test: [%d x %d] samples and genes.\n', size(dataset_info.DatasetTr.Gene_Expression), size(dataset_info.DatasetTe.Gene_Expression))
    
    %% Loop over methods
    for mi=1:n_meth
        fprintf(['--- Evaluate [%s, %s]' repmat('-',1,40) '\n'], net_lst{ni}, method_lst{mi});
        result_name = sprintf([result_path '%s_MSN-%03d_MTN-%s.mat'], dataset_list(1).name(1:end-4), opt_info.MAX_N_SUBNET, method_lst{mi});
        if exist(result_name, 'file')
            fprintf('[i] Results are already computed --> [%s] \n', result_name);
            result = load(result_name);
            te_auc(ni, mi) = result.te_auc;
            fprintf('AUC was [%0.2f]\n', te_auc(ni, mi)*100);
            continue;
        end
        
        %% Evaluate model
        fprintf('[%d/%d] Evaluating [%s] model on [%s] ...\n', mi, n_meth, method_lst{mi}, datestr(now));
        switch method_lst{mi}
            case 'iPark'
                result = perf_iPark(dataset_info, opt_info, 'Mean');
            case 'RI-iPark'
                result = perf_iPark(dataset_info, opt_info, 'RI');
            case 'DA2Lex'
                InfSN_info = opt_info;
                InfSN_info.MAX_N_SUBNET = inf;
                result = perf_DA2Lex(dataset_info, InfSN_info, 'Mean');
            case {'GLasso2' 'GLasso3' 'GLasso5' 'GLasso7' 'GLasso10' 'GLasso15' 'GLasso20'}
                opt_gls = opt_info;
                opt_gls.lam_list = [zeros(20,1) logspace(log10(1e-2), 0, 20)'];
                opt_gls.MAX_SUBNET_SIZE = str2double(method_lst{mi}(7:end));
                result = perf_GLasso(dataset_info, opt_gls);
            case 'NetGL'
                opt_ngl = opt_info;
                opt_ngl.lam_list = [zeros(20,1) logspace(log10(1e-2), 0, 20)'];
                tmp_res_ptr = sprintf('./Results_Files/DID_%s_*_MSN-500_MTN-NetLasso.mat', ds_id);
                tmp_res_info = dir(tmp_res_ptr);
                if numel(tmp_res_info)~=1, error(); end
                fprintf('[i] Net lasso result is found in [%s], loading ...\n', tmp_res_info.name);
                tmp_info = load(sprintf('./Results_Files/%s', tmp_res_info.name));
                opt_ngl.MAX_N_Gene = tmp_info.BestNetwork;
                clear tmp_info tmp_res_info tmp_res_ptr
                result = perf_NetGL(dataset_info, opt_ngl);
            case 'CFGLasso'
                result = perf_CFGLasso(dataset_info, opt_info);
            case {'FERALAvg' 'FERALAvgStdInt' 'FERALAvgStd' 'FERALInt'}
                feral_info = opt_info;
                feral_info.CompFeat = arrayfun(@(i) method_lst{mi}(i:i+2), 6:3:numel(method_lst{mi}), 'UniformOutput', 0);
                result = perf_FERAL(dataset_info, feral_info);
            case {'Lasso'}
                result = perf_Lasso(dataset_info, opt_info);
            case {'NetLasso'}
                result = perf_NetLasso(dataset_info, opt_info);
            case 'LExAG'
                result = perf_LExAG(dataset_info, opt_info);
            case 'TLEx'
                result = perf_TLEx(dataset_info, opt_info);
            case 'TAgLEx'
                result = perf_TAgLEx(dataset_info, opt_info);
            case 'TRgLEx'
                opt_rndg = opt_info;
                opt_rndg.UseRndGene = 1;
                result = perf_TAgLEx(dataset_info, opt_rndg);
            case {'NMC' 'TNMC'}
                opt_nmc = opt_info;
                if strcmpi(method_lst{mi}(1), 'T')
                    opt_nmc.UseTTest = 1;
                end
                result = perf_TTNMC(dataset_info, opt_nmc);
            case 'TNMCAd'
                result = perf_TTNMC(dataset_info, setfield(opt_info, 'FindK', 1));
            case 'TAgNMC'
                result = perf_TAgNMC(dataset_info, opt_info);
            case 'TRgNMC'
                opt_rndg = opt_info;
                opt_rndg.UseRndGene = 1;
                result = perf_TAgNMC(dataset_info, opt_rndg);
            case 'Regress'
                result = perf_Regress(dataset_info, opt_info);
            case 'TReg'
                result = perf_Regress(dataset_info, struct('K', MAX_N_SUBNET));
            case 'RegAG'
                result = perf_RegAG(dataset_info, opt_info);
            case {'KNN0','KNN1','KNN3','KNN5','KNN7','TKNN0'}
                opt_knn = opt_info;
                opt_knn.K = str2double(method_lst{mi}(end));
                if strcmpi(method_lst{mi}(1), 'T')
                    opt_knn.UseTTest = 1;
                end
                result = perf_KNN(dataset_info, opt_knn);
            case {'NB' 'TNB'}
                opt_nb = opt_info;
                if strcmpi(method_lst{mi}(1), 'T')
                    opt_nb.UseTTest = 1;
                end
                result = perf_NB(dataset_info, opt_nb);
            case {'LDA' 'TLDA'}
                opt_lda = opt_info;
                if strcmpi(method_lst{mi}(1), 'T')
                    opt_lda.UseTTest = 1;
                end
                result = perf_LDA(dataset_info, opt_lda);
            case {'NN' 'TNN'}
                opt_nn = opt_info;
                if strcmpi(method_lst{mi}(1), 'T')
                    opt_nn.UseTTest = 1;
                end
                opt_nn.GridSearch = 1;
                result = perf_NN(dataset_info, opt_nn);
            case {'SVM-Lin','SVM-RBF','TSVM-Lin','TSVM-RBF'}
                opt_svm = opt_info;
                if strcmpi(method_lst{mi}(end-2:end), 'rbf')
                    opt_svm.kernel_name = 'rbf';
                else
                    opt_svm.kernel_name = 'linear';
                end
                if strcmpi(method_lst{mi}(1), 'T')
                    opt_svm.UseTTest = 1;
                end
                opt_svm.GridSearch = 1;
                result = perf_SVM(dataset_info, opt_svm);
            case {'DT' 'TDT'}
                opt_dt = opt_info;
                if strcmpi(method_lst{mi}(1), 'T')
                    opt_dt.UseTTest = 1;
                end
                opt_dt.GridSearch = 1;
                result = perf_DecisionTree(dataset_info, opt_dt);
            case {'RnFrst' 'TRnFrst'}
                opt_rf = opt_info;
                if strcmpi(method_lst{mi}(1), 'T')
                    opt_rf.UseTTest = 1;
                end
                opt_rf.GridSearch = 1;
                result = perf_RandomForest(dataset_info, opt_rf);
            otherwise
                error('Unknown Method.');
        end
        te_auc(ni, mi) = result.te_auc;
        result.Method_Name = method_lst{mi};
        result.Dataset_Name = dataset_name;
        %if ismac, return; end
        
        %% Saving results
        fprintf('Saving results in [%s].\n', result_name);
        save(result_name, '-struct', 'result');
        clear result
        fprintf('\n');
    end
end

%% Print results
te_auc(end+1, :) = mean(te_auc, 1, 'omitnan');
te_auc(:, end+1) = mean(te_auc, 2, 'omitnan');
fprintf([repmat('$',1,50) '\n']);
fprintf('Final results:\n');
disp(array2table(te_auc, 'RowNames', [strrep(net_lst,'-','_') 'AVG'] , 'VariableNames', [strrep(method_lst,'-','_') 'AVG']));
fprintf([repmat('$',1,50) '\n']);

