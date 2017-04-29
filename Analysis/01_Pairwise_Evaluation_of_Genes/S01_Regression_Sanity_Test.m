
n_sample = 5000;
data = randn(n_sample, 5);
r = [1 1.1 0.5 -1 0.5]'; % Only two nonzero coefficients
label = ((data*r + randn(n_sample,1)*7)>0)*2-1;

data = data + randn(n_sample,5) * 2;
iTr = false(n_sample,1); iTr(randperm(n_sample, n_sample*0.7)) = 1;
iTe = ~iTr;
zTr = data(iTr, :);
zTe = data(iTe, :);
lTr = label(iTr, :);
lTe = label(iTe,:);

cvp = cvpartition(lTr, 'KFold', 4);
for fi=1:cvp.NumTestSets
	Fold_Index(cvp.test(fi),1) = fi;
end

[las_B, las_fit] = lasso(zTr, lTr, 'CV', cvp, 'NumLambda', 20);

lasso_opt = {'lassoType', 't', 'CV', 1, 'iCvPar', Fold_Index, 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1};
[me_B, me_fit] = lassoEx(zTr, lTr, lasso_opt{:});

imagesc([zscore(las_B); zscore(me_B)])

figure();
plot(las_fit.MSE', 1 - me_fit.MSE', '*');