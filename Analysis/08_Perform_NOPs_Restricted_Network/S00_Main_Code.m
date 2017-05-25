function S00_Main_Code(Target_Study, Target_Repeat, method_lst, net_lst, MAX_N_SUBNET)
%% Run 
%{
for ri in `seq 1 10`; do
for si in `seq 1 12`; do
PARAM="$si,$ri"; sbatch --job-name=NE-$PARAM --output=Logs/NE-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=8GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S00_Main_Code "$PARAM";
done;
read -p "Press a key" -t 1800
done

UMC: PARAM="si,$ri,{'TAgNMC','TNMC','TLEx','TAgLEx'},{'Random-T00010'},10"; qsub -N "NE-$PARAM" ~/bulk/env/run_Matlab.sh S00_Main_Code "$PARAM";
%}

%% ####
if ismac
	fprintf('*** Warning!: Running on debug mode.\n');
	Target_Repeat = 1;
	Target_Study = 1;
	method_lst = {'TAgNMC'};
	net_lst = {'Random-T00050'};
	MAX_N_SUBNET = 20;
end

%% Initialization
clc
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath('../_Utilities/');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
dataset_path = './Dataset_Files/';
result_path = './Results_Files/';
opt_info.lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 0};
if ~exist('MAX_N_SUBNET', 'var')
	opt_info.MAX_N_SUBNET = 500;
else
	opt_info.MAX_N_SUBNET = MAX_N_SUBNET;
end

if ~exist('net_lst', 'var') || isempty(net_lst)
	net_lst = {
		'CrSyn-T10000', 'CrMinSyn-T10000' ... % 'SyNet-T10000', 'MinSyn-T10000', 'AvgSyn-T10000', 
		'AbsCorr-T10000','Random-T10000', ... % 'Corr-T10000',
		'STRING-T10000', 'KEGG-T10000'
		};
end
if ~exist('method_lst', 'var') || isempty(method_lst)
	method_lst = {'TLEx' 'Lasso' 'Feral' 'RI-Feral' 'CFGLasso' 'iPark' 'RI-iPark' 'LExAG'}; % 'Regress' 'AS-Feral'
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
	ds_id = S02_GenerateDataset(cv_id, net_lst{ni});
	
	%% Load Dataset info
	dataset_list = dir([dataset_path 'DID_' ds_id '*.mat']);
	if numel(dataset_list)~=1, error('Missing or duplicate dataset found.. \n[%s]\n', strjoin({dataset_list.name}, ', ')); end
	dataset_name = [dataset_path dataset_list(1).name];
	fprintf('Loading dataset [%s] ...\n', dataset_name);
	dataset_info = load(dataset_name);
	
	%% Loop over methods
	for mi=1:n_meth
		fprintf(['--- Evaluate [%s, %s]' repmat('-',1,40) '\n'], net_lst{ni}, method_lst{mi});
		result_name = sprintf([result_path '%s_MSN-%03d_MTN-%s.mat'], dataset_list(1).name(1:end-4), opt_info.MAX_N_SUBNET, method_lst{mi});
		if exist(result_name, 'file')
			fprintf('[i] Results are already computed --> [%s] \n', result_name);
			continue;
		end
		
		%% Evaluate model
		fprintf('[%d/%d] Evaluating [%s] model ...\n', mi, n_meth, method_lst{mi});
		switch method_lst{mi}
			case 'iPark'
				result = perf_iPark(dataset_info, opt_info, 'Mean');
			case 'RI-iPark'
				result = perf_iPark(dataset_info, opt_info, 'RI');
			case 'iChuang'
				result = perf_iChuang(dataset_info, opt_info);
			case 'iTaylor'
				result = perf_iTaylor(dataset_info, opt_info);
			case 'Feral'
				result = perf_Feral(dataset_info, opt_info, 'Mean');
			case 'RI-Feral'
				result = perf_Feral(dataset_info, opt_info, 'RI');
			case 'AS-Feral'
				highSN_info = opt_info;
				highSN_info.MAX_N_SUBNET = 10000;
				result = perf_Feral(dataset_info, highSN_info, 'RI');
			case 'GLasso'
				result = perf_GLasso(dataset_info, opt_info);
			case 'CFGLasso'
				result = perf_CFGLasso(dataset_info, opt_info);
			case 'Lasso'
				result = perf_Lasso(dataset_info, opt_info);
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
			case 'TNMC'
				result = perf_TTNMC(dataset_info, opt_info);
			case 'TAgNMC'
				result = perf_TAgNMC(dataset_info, opt_info);
			case 'TRgNMC'
				opt_rndg = opt_info;
				opt_rndg.UseRndGene = 1;
				result = perf_TAgNMC(dataset_info, opt_rndg);
			case 'Regress'
				result = perf_Regress(dataset_info, opt_info);
			case 'RegAG'
				result = perf_RegAG(dataset_info, opt_info);
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

