
%% Initialization
addpath('../../../../Useful_Sample_Codes/ErrorbarEx/');
GEPath = getPath('SyNet');
opt_info.lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 0};

%% Load data
load(GEPath, 'Gene_Expression', 'Patient_Label', 'Gene_Name', 'Study_Index');
zData = zscore(Gene_Expression);
Patient_Label = double(Patient_Label==1)*2-1;
n_study = max(Study_Index);
n_gene = size(zData, 2);

%% Training
cr_mat = zeros(n_study, n_gene);
for si=1:n_study
	iTr = Study_Index == si;
	cr_mat(si,:) = corr(zData(iTr,:), Patient_Label(iTr), 'Type', 'Spearman');
end

%% Sort
length_lst = [10 20 50 100];
n_len = numel(length_lst);
auc_mat = zeros(n_study, n_len);
for li=1:n_len
	for si=1:n_study
		n_feat = length_lst(li);
		rest_ind = setdiff(1:n_study, si);
		avg_cr = mean(cr_mat(rest_ind,:),1);
		std_cr =  std(cr_mat(rest_ind,:), 0, 1);
		std_cr(abs(avg_cr)<0.05) = nan;
		[val, ind] = sort(std_cr, 'Ascend');
		gind = ind(1:n_feat);
		
		%{
			plot(1:n_feat, avg_cr(gind), 'b');
			hold on
			errorbarEx(1:n_feat, avg_cr(gind), std_cr(gind), std_cr(gind), 2, 0.2, [1 0 0]);
		%}
		iTr = Study_Index ~= si;
		iTe = Study_Index == si;
		zTr = zData(iTr, gind);
		zTe = zData(iTe, gind);
		lTr = Patient_Label(iTr);
		lTe = Patient_Label(iTe);
		result = LassoWithCV(@lassoEx, zTr, lTr, zTe, lTe, Study_Index(iTr), opt_info.lasso_opt);
		auc_mat(si, li) = result.te_auc;
		fprintf('@@@@@ test performance [%0.2f%%] AUC.\n', result.te_auc*100);
	end
end
mean(auc_mat)
plot(auc_mat);
legend(cellstr(num2str(length_lst')))
hold on
load('../01_Pairwise_Evaluation_of_Genes/Baseline_AUCs/BA_CV01_TAgNMC_Random-T00020.mat');
plot(BaseLine_AUC, 'LineWidth', 2);