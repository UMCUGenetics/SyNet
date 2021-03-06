function S00_Main_Code(Target_Study, Target_Repeat, method_lst, net_lst, MAX_N_SUBNET)
%% Run 
%{
for ri in `seq 1 10`; do
for si in `seq 1 14`; do
PARAM="$si,$ri"; sbatch --job-name=NE-$PARAM --output=Logs/NE-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=10GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S00_Main_Code "$PARAM";
done;
read -p "Press a key" -t 1800
done

UMC: PARAM="$si,$ri,{'TAgNMC','TNMC','TLEx','TAgLEx'},{'Random-T00010'},10"; qsub -N "NE-$PARAM" ~/bulk/env/run_Matlab.sh S00_Main_Code "$PARAM";
%}

%% ####
if ismac || ispc
	fprintf('*** Warning!: Running on debug mode.\n');
	Target_Study = 3;
    Target_Repeat = 1;
	method_lst = {'TReg'};
	net_lst = {'None-G11748'};
	MAX_N_SUBNET = 500;
end

%% Initialization
clc
addpath(genpath('../../../../Useful_Sample_Codes/ShowProgress'));
addpath('../_Utilities/');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
addpath(genpath('../../../../Useful_Sample_Codes/getAUC'));
dataset_path = './Dataset_Files/';
result_path = './Results_Files/';
opt_info.lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};
if ~exist('MAX_N_SUBNET', 'var')
	opt_info.MAX_N_SUBNET = 500;
else
	opt_info.MAX_N_SUBNET = MAX_N_SUBNET;
end

if ~exist('net_lst', 'var') || isempty(net_lst)
	net_lst = {
        'None-G11748'
		'AvgSynACr-P10000','AvgSyn-P10000', ...
		'AbsCorr-P10000', ...
		'STRING-P10000','KEGG-P10000','Random-P10000','I2D-P10000','HPRD-P10000','MSigDB-P10000', ...
		'AvgSynACr-G00500','AvgSyn-G00500', ...
		'AbsCorr-G00500', ...
		'STRING-G00500','KEGG-G00500','Random-G00500','I2D-G00500','HPRD-G00500','MSigDB-G00500', ...
		};
end
if ~exist('method_lst', 'var') || isempty(method_lst)
	method_lst = {'TNMC' 'TNMCAd' 'TLEx' 'Lasso' 'GLasso' 'CFGLasso'}; % 'Regress' 'AS-Feral'
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
	fprintf('Dataset has Train: [%d,%d], Test: [%d,%d] samples and genes.\n', size(dataset_info.DatasetTr.Gene_Expression), size(dataset_info.DatasetTe.Gene_Expression))
	
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
			case {'TNMC0' 'TNMC20' 'TNMC50' 'TNMC100' 'TNMC500'}
				result = perf_TNMC(dataset_info, struct('K', str2double(method_lst{mi}(5:end))));
            case {'KNN0','KNN1','KNN3','KNN5','KNN7'}
                result = perf_KNN(dataset_info, setfield(opt_info, 'K', str2double(method_lst{mi}(4:end))));
			case 'Lasso'
				result = perf_Lasso(dataset_info, opt_info);
			case 'LExAG'
				result = perf_LExAG(dataset_info, opt_info);
			case 'TLEx'
				result = perf_TLEx(dataset_info, opt_info);
			otherwise
				error('Unknown Method.');
		end
		te_auc(ni, mi) = result.te_auc;
		result.Method_Name = method_lst{mi};
		result.Dataset_Name = dataset_name;
		
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

