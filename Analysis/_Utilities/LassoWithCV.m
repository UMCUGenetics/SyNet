function result = LassoWithCV(func, xTr, lTr, xTe, lTe, Study_Index, lasso_opt)

%% Initialization
if ~exist('lasso_opt', 'var')
	lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-2, 'n_lC', 20, 'lC_ratio', 1e-2, 'verbose', 1};
end
if ~exist('func', 'var') || isempty(func)
	func = @lassoEx;
end
[n_Tr, n_gene] = size(xTr);

%% Traning
fprintf('Training [%s] over all samples [%d x %d]...\n', func2str(func), n_Tr, n_gene);
[opt_B, opt_fit] = func(xTr, lTr, lasso_opt{:});

%% Showing results
n_lam = size(opt_B, 2);
tr_auc_lam = zeros(1, n_lam);
te_auc_lam = zeros(1, n_lam);
for i=1:n_lam
	tr_auc_lam(i) = getAUC(lTr, xTr*opt_B(:,i), 50);
	if ~isempty(xTe)
		te_auc_lam(i) = getAUC(lTe, xTe*opt_B(:,i), 50);
	end
end
fprintf('    Index: '); fprintf('%5d  ', 1:n_lam); fprintf('\n');
fprintf('Train AUC: '); fprintf('%0.3f, ', tr_auc_lam); fprintf('\n');
fprintf(' Test AUC: '); fprintf('%0.3f, ', te_auc_lam); fprintf('\n');

%% Preparing CV-Index
[Study_lst, ~, Fold_Index] = unique(Study_Index, 'Stable');
if numel(Study_lst)==1
	fprintf('[i] Warning: Only one study is found in the training set. Running a 5-CV instead.\n');
	cv_obj = cvpartition(lTr, 'kFold', 5);
	for fi=1:cv_obj.NumTestSets
		Fold_Index(cv_obj.test(fi), 1) = fi;
	end
end
n_iFold = max(Fold_Index);

%% Cross-validation
fold_auc = zeros(n_iFold, n_lam);
fprintf('Training over [%d] folds ...\n', n_iFold);
for fi=1:n_iFold
	fTr = Fold_Index~=fi;
	fTe = Fold_Index==fi;
	fprintf('Fold [%02d/%02d]: Training over [#Tr=%d, #Te=%d], ', fi, n_iFold, sum(fTr), sum(fTe));
	[fold_B, ~] = func(xTr(fTr,:), lTr(fTr), lasso_opt{:});
	if size(fold_B,2)<n_lam
		fold_B = [fold_B zeros(n_gene, n_lam-size(fold_B,2))];
	end
	fold_pred = xTr(fTe,:)*fold_B;
	for li=1:n_lam
		fold_auc(fi,li) = getAUC(lTr(fTe), fold_pred(:,li), 50);
	end
	[~, fold_MIN_MSE_IND] = min(1 - fold_auc(fi,:));
	fprintf('Smallest MSE is on [%d]\n', fold_MIN_MSE_IND);
end
opt_fit.MSE = 1 - mean(fold_auc, 1);
[opt_fit.MinMSE, opt_fit.IndexMinMSE] = min(opt_fit.MSE);
fprintf('Best Lambda is determined as [%d].\n', opt_fit.IndexMinMSE);
fprintf('Number of non-zero elements is: [%d]\n', sum(abs(opt_B(:,opt_fit.IndexMinMSE))>0));

%% Saving
result.TrLbl = lTr;
result.TeLbl = lTe;
result.predTrLbl = xTr * opt_B(:, opt_fit.IndexMinMSE);
result.predTeLbl = xTe * opt_B(:, opt_fit.IndexMinMSE);
result.B = opt_B;
result.fit = opt_fit;
result.tr_auc_lam = tr_auc_lam;
result.te_auc_lam = te_auc_lam;
result.tr_auc = tr_auc_lam(opt_fit.IndexMinMSE);
result.te_auc = te_auc_lam(opt_fit.IndexMinMSE);
end