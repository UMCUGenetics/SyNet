function result = Perf_NMC(Dataset_info)

%% Initialization

%% Prepare data
zTr = Dataset_info.zTr;
zTe = Dataset_info.zTe;
lTr = Dataset_info.lTr;
lTe = Dataset_info.lTe;

%% Train NMC
fprintf('Training NMC ...\n');
pred = nmc(zTr, lTr, zTe);

%% Evaluation
fprintf('Evaluation ...\n');
result.te_auc = getAUC(lTe, pred, 50);
end