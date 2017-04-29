
rng default % For reproducibility
X = randn(100,5);
r = [0;2;0;-3;0]; % Only two nonzero coefficients
Y = X*r + randn(100,1)*5; % Small added noise
l = (Y>0.5)*2-1;
[las_B, las_fit] = lasso(X,l, 'NumLambda', 20);

lasso_opt = {'lassoType', 't', 'CV', [], 'relTol', 5e-1, 'n_lC', 20, 'lC_ratio', 1e-2, 'paroptions', statset('UseParallel',false), 'verbose', 1};
[ex_B, ex_fit] = lassoEx(X+randn(100,5),l, lasso_opt{:});

imagesc([las_B; ex_B]);
colormap(jet);
colorbar

