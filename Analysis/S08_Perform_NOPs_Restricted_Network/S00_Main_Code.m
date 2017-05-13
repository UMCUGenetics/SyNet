function S00_Main_Code(Target_Study, Target_Repeat, method_lst, net_lst)
%% Run 
%{
for ri in `seq 1 10`; do
for si in `seq 1 12`; do
PARAM="$si,$ri"; sbatch --job-name=NE-$PARAM --output=Logs/NE-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=8GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S00_Main_Code "$PARAM";
done;
read -p "Press a key" -t 1800
done
%}

%% ####
if ismac
	fprintf('*** Warning!: Running on debug mode.\n');
	Target_Repeat = 1;
	Target_Study = 1;
	net_lst = {'CrSyn-T10000'};
	method_lst = {'RI-Feral'};
end

%% Initialization
clc
cv_id = sprintf('CVT01_Si%02d-Ri%03d', Target_Study, Target_Repeat);
if ~exist('net_lst', 'var') || isempty(net_lst)
	net_lst = {
		'SyNet-T10000' ...
		'Corr-T10000','AbsCorr-T10000','Random-T10000', ...
		};
end
if ~exist('method_lst', 'var') || isempty(method_lst)
	method_lst = {'Regress' 'AS-Feral' 'TLEx' 'Feral'  'RI-Feral' 'iPark' 'RI-iPark' 'Lasso' 'LExAG'}; 
end
n_net = numel(net_lst);
n_meth = numel(method_lst);

%% Main Loop
fprintf([repmat('/',1,20) ' Start of main loop ' repmat('/',1,20) '\n']);
fprintf('[i] CV ID is: %s\n', cv_id);
fprintf('[i] Method list is: %s\n', strjoin(method_lst, ', '));
fprintf('[i] Network list is: %s\n', strjoin(net_lst, ', '));
fprintf([repmat('/',1,60) '\n']);

te_auc = nan(n_net, n_meth);
for ni=1:n_net %% Dataset = Data + Network
	fprintf(['****** Network Gen [%s] ' repmat('*',1,40) '\n'], net_lst{ni});
	ds_id = S02_GenerateDataset(cv_id, net_lst{ni});
	
	for mi=1:n_meth %% Methods
		fprintf(['--- Evaluate [%s, %s]' repmat('-',1,40) '\n'], net_lst{ni}, method_lst{mi});
		[~, te_auc(ni, mi)] = S03_Evaluating_Models(method_lst(mi), ds_id);
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

