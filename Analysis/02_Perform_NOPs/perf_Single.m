function result = perf_Single(dataset_info, opt_info)

%% Initialization
xTr = dataset_info.DatasetTr.Gene_Expression;
lTr = dataset_info.DatasetTr.Patient_Label;
xTe = dataset_info.DatasetTe.Gene_Expression;
lTe = dataset_info.DatasetTe.Patient_Label;
n_gene = size(xTr, 1);
MAX_N_SUBNET = opt_info.MAX_N_SUBNET;

%% Normalization
fprintf('Normalizing data ...\n');
zTr = zscore(xTr);
zTe = zscore(xTe);

%% Traning the final model
fprintf('Training the final model over [%d] features...\n', n_gene);
[opt_B, opt_fit] = lassoEx(zTr, lTr, opt_info.lasso_opt{:}, 'iCvPar', dataset_info.DatasetTr.iCvPar);
fprintf('Final training is done. [%d] non-zero features identified.\n', sum(abs(opt_B(:, opt_fit.IndexMinMSE))>0));

%% Evaluating the model
vec_B = opt_B(:, opt_fit.IndexMinMSE);
[~, topB_ind] = sort(abs(vec_B), 'Descend');
vec_B(topB_ind(MAX_N_SUBNET+1:end)) = 0;

g = zTr*vec_B;
tr_auc = getAUC(lTr, g, 50);

g = zTe*vec_B;
te_auc = getAUC(lTe, g, 50);

fprintf('@@@@@ Final test performance for this dataset is [%0.2f%%] AUC.\n', te_auc*100);

result.B = opt_B;
result.fit = opt_fit;
result.tr_auc = tr_auc;
result.te_auc = te_auc;

% for i=1:size(zTr,2) %###
% 	tr_auc(i,1) = getAUC(lTr, zTr(:,i), 50);
% 	te_auc(i,1) = getAUC(lTe, zTe(:,i), 50);
% end
% plot(tr_auc, 'b'); hold on
% plot(te_auc, 'r');
end


