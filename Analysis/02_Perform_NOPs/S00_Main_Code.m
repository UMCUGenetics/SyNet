function S00_Main_Code(data_lst, method_lst, net_lst, cv_id)
%{
for fn in CV_Files/*.mat; do
	fi=${fn:13:14}
	echo -e "///\nFile: $fn"
	echo "ID: $fi"
	si=${fi:8:2}
	sname=SyNet
	fparam="{},[],{'Corr','AbsCorr','Random','String','Multinet','KEGG','MSigDB','DSN-$sname$si-P99.90','DSN-$sname$si-P99.99','DSN-$sname$si-P99.999','DSN-$sname$si-T00500','DSN-$sname$si-T01000','DSN-$sname$si-T05000','DSN-$sname$si-T10000','DSN-$sname$si-T20000'},'$fi'"; 
	sbatch --job-name=$fi --output=Logs/Main-$fi.%J_%a-%N.out --partition=general --qos=long --mem=16GB --time=7-00:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S00_Main_Code $fparam
done
%}
	
clc;
% if ismac
% 	fprintf('*** Warning!: Running on debug mode.\n');
% 	method_lst = {'iChuang'};
% 	data_lst = {};
% 	net_lst = {'DSN-SyNetS1-T500'};
% 	cv_id = '170411001004';
% end

%% Initialization
if isempty(net_lst)
	net_lst = {
		'Corr','AbsCorr','Random','String','Multinet','KEGG','MSigDB',...
		'DSN-SyNet-P99.90','DSN-SyNet-T10000' ...
		};
end
n_net = numel(net_lst);

if isempty(method_lst)
	method_lst = {'iPark', 'iChuang', 'iTaylor', 'Feral'}; 
end
n_meth = numel(method_lst);

if isempty(data_lst)
	if isempty(cv_id)
		error('CV ID is required when data names are not given.\n'); 
	end
	fprintf('[i] CV ID is defined to be [%s]\n', cv_id);
	cv_path = './CV_Files/';
	cv_info = dir([cv_path 'CID-' cv_id '*.mat']);
	cv_name = [cv_path cv_info.name];
	fprintf('CV file found [%s]\n', cv_name);
	data_lst = regexp(cv_name, '.*_ETr-(\w+)_ETe-(\w+).*', 'tokens');
	data_lst = {data_lst{1}{1} data_lst{1}{2}};
	fprintf('Data list is now set on: %s\n', strjoin(data_lst, ', '));
end

% if any(~cellfun(@isempty, strfind(net_lst, data_lst{1}))) || ...
%    any(~cellfun(@isempty, strfind(net_lst, data_lst{2})))
% 	error('Trained networks can not be used as test or training set.');
% end

%% CV Generation
% if isempty(cv_id)
% 	cv_id = S01_GenerateCVSet(data_lst{1}, data_lst{2});
% 	fprintf('[i] CV ID is set to [%s]\n', cv_id);
% end

%% Main Loop
fprintf([repmat('/',1,20) ' Start of main loop ' repmat('/',1,20) '\n']);
fprintf('[i] Data list is: %s\n', strjoin(data_lst, ', '));
fprintf('[i] Method list is: %s\n', strjoin(method_lst, ', '));
fprintf('[i] Network list is: %s\n', strjoin(net_lst, ', '));
fprintf([repmat('/',1,60) '\n']);

te_auc = nan(n_net, n_meth);
for ni=1:n_net %% Dataset = Data + Network
	fprintf(['****** Network Gen [%s] ' repmat('*',1,40) '\n'], net_lst{ni});
	ds_id = S02_GenerateDataset(cv_id, net_lst{ni});
	
	for mi=1:n_meth %% Methods
		fprintf(['--- Evaluate [%s, %s, %s, %s]' repmat('-',1,40) '\n'], data_lst{1}, data_lst{2}, net_lst{ni}, method_lst{mi});
		[~, te_auc(ni, mi)] = S03_Evaluating_Models(method_lst(mi), ds_id);
		fprintf('\n');
	end
end

%% Print results
te_auc(end+1, :) = mean(te_auc, 1, 'omitnan');
te_auc(:, end+1) = mean(te_auc, 2, 'omitnan');
fprintf([repmat('$',1,50) '\n']);
fprintf('Training GE: %s,\tTesting GE: %s\n\n', data_lst{1}, data_lst{2});
disp(array2table(te_auc, 'RowNames', [net_lst 'AVG'] , 'VariableNames', [method_lst 'AVG']));
fprintf([repmat('$',1,50) '\n']);
