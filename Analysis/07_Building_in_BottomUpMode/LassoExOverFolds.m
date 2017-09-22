function result = LassoExOverFolds(zTr, lTr, zTe, lTe, Study_Index)

%% Prepare data
n_gene = size(zTr, 2);
[~, ~, Fold_Index] = unique(Study_Index, 'Stable');
n_iFold = max(Fold_Index);

%% Combination
B_fold = zeros(n_gene, n_iFold);
for fi=1:n_iFold
	fTe = Fold_Index==fi;
	fprintf('Fold [%02d/%02d]: Training over [%d] samples ...\n', fi, n_iFold, sum(fTe));
	B_fold(:, fi) = regress(lTr(fTe), zTr(fTe,:));
end
opt_B = mean(B_fold, 2);

%% Evaluation
tr_auc = getAUC(lTr, zTr*opt_B, 50);
te_auc = getAUC(lTe, zTe*opt_B, 50);
fprintf('AUC is: Tr=%0.3f, Tr=%0.3f\n', tr_auc, te_auc);

%% Saving
result.B = opt_B;
result.tr_auc = tr_auc;
result.te_auc = te_auc;
end