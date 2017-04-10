function S03_PerformPWR(ge_name, batch_be, batch_en)
% Run: S03_PerformPWR('SyNet', 10, 20)

%% Initialization
addpath('../../../../Useful_Sample_Codes/getAUC/');
addpath('../../../../Useful_Sample_Codes/fastAUC/');
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
pwr_path = ['./PWR_Files/' ge_name '/'];
cv_path = './CV_Files/';
[~, ~] = mkdir(pwr_path);
fprintf('/// Running pairwise evaluation for [%s], we will test [%d-%d] pairs.\n', ge_name, batch_be, batch_en);

%% Load data
GeneExpression_Path = getPath(ge_name);
fprintf('Loading [%s]\n', GeneExpression_Path);
data = load(GeneExpression_Path, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
zData = zscore(data.Gene_Expression);
Patient_Label = double(data.Patient_Label);
Gene_Name = data.Gene_Name;
[n_sample, n_gene] = size(zData);
fprintf('Data loaded [n_sample=%d, n_gene=%d]...\n', n_sample, n_gene);
clear data

%% Identify pairs
fprintf('Loading pair list.\n');
[gi, gj] = find(triu(ones(n_gene), 0));
pair_list = [gi gj];
clear gi gj
n_total = size(pair_list,1);
fprintf('In total [%d] gene pairs exist.\n', n_total);
if batch_en > n_total, batch_en = n_total; end
pair_list = pair_list(batch_be:batch_en, :);
n_pair = size(pair_list, 1);
fprintf('Running pairwise check for [%d] pairs.\n', n_pair);

%% Load CV info
cv_name = [cv_path 'CV_' ge_name '_CVT01.mat'];
cv_info = load(cv_name);
if ~isequal(cv_info.Patient_Label, Patient_Label), error(); end
cv_obj = cv_info.cv_obj;
[n_fold, n_rep] = size(cv_obj);
fprintf('Loading CV info from [%s], [%d] folds and [%d] repeats are loaded.\n', cv_name, n_fold, n_rep);

%% Main loop
fprintf('Pairwise comparison started at: %s\n', datetime);
fprintf('Checking pairs:\n');
auc_pair = zeros(n_pair, 2+n_fold);
for pi=1:n_pair
	showprogress(pi, n_pair);
	if mod(pi, 1000)==0
		fprintf('Currently at [%08d/%08d] for genes [%05d/%05d]\n', pi, n_pair, pair_list(pi,:));
	end
	auc_rep = zeros(n_fold, n_rep);
	for ri=1:n_rep
		for fi=1:n_fold
			%fprintf('Rep [%d], fold [%d]\n', ri, fi);
			iTr = cv_obj(fi, ri).iTr;
			iTe = cv_obj(fi, ri).iTe;
			zTr = zData(iTr, pair_list(pi,:));
			lTr = Patient_Label(iTr);
			zTe = zData(iTe, pair_list(pi,:));
			lTe = Patient_Label(iTe);
			
			if pair_list(pi,1)==pair_list(pi,2)
				pred = zTe(:,1);
			else
				B = regress(lTr, zTr);
				pred = zTe * B;
			end
			% auc = getAUC(lTe, pred);
			% [~,~,~,auc] = perfcurve(lTe, pred, 1)
			tmp = fastAUC(lTe, pred, 1); auc_rep(fi, ri) = max([1-tmp tmp]);
		end
	end
	auc_pair(pi, :) = [pair_list(pi,:) mean(auc_rep, 2)'];
end
fprintf('Pairwise comparison finished at: %s\n', datetime);

%% Saving
sav_name = sprintf('%sPWR_%s_%08d-%08d.mat', pwr_path, ge_name, batch_be, batch_en);
fprintf('Saving result in [%s]', sav_name);
save(sav_name, 'auc_pair', 'Patient_Label', 'cv_name', 'cv_obj', 'pair_list', 'Gene_Name');
fprintf('Process finished at: %s\n', datetime);
end

