function S04_Evaluate_Integration_Methods(Target_Study, Target_Repeat, method_name, net_name)

%% Initialization
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
addpath('../../../../Useful_Sample_Codes/getAUC/');
addpath('../_Utilities/');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
tr_name = 'SyNet';
te_name = 'SyNet';
cv_ind = 50;
if ispc
    method_name = 'PCA1';
    net_name = 'Random-NN05';
    Target_Study = 3;
    Target_Repeat = 2;
end

%% Load indices
cv_name = sprintf('./CV_Files/CV_%s-%s_CVT%02d_Si%02d-Ri%03d.mat', tr_name, te_name, cv_ind, Target_Study, Target_Repeat);
fprintf('Loading the CV file: %s\n', cv_name);
cv_info = load(cv_name);
cv_info.cv_name = cv_name;

%% Load Train data
fprintf('Loading Train expression data from: %s\n', cv_info.tr_info.GEPath);
data = load(cv_info.tr_info.GEPath, 'Gene_Expression', 'Patient_Label', 'Study_Index', 'Gene_Name');
if ~isequal(data.Gene_Name, cv_info.tr_info.Gene_Name), error(); end
iTr = cv_info.tr_info.CVInd;
zTr = zscore(data.Gene_Expression(iTr,:));
lTr = (data.Patient_Label(iTr)==1)*2-1;

%% Load Test data
fprintf('Loading Test expression data from: %s\n', cv_info.te_info.GEPath);
data = load(cv_info.te_info.GEPath, 'Gene_Expression', 'Patient_Label', 'Study_Index', 'Gene_Name');
if ~isequal(data.Gene_Name, cv_info.te_info.Gene_Name), error(); end
iTe = cv_info.te_info.CVInd;
zTe = zscore(data.Gene_Expression(iTe,:));
lTe = (data.Patient_Label(iTe)==1)*2-1;
clear data

%% Load neighbors
nei_name = ['./NetNei_Files/NetNei_' net_name '.mat'];
fprintf('Loading neighbor file: %s\n', nei_name);
nei_info = load(nei_name);
if ~isequal(nei_info.Gene_Name, cv_info.tr_info.Gene_Name), error(); end
n_epoch = numel(nei_info.SubNet_Full);
nei_info.nei_name = nei_name;

%% Main loop
fprintf('Comparing performance for [%d] neighbor sets:\n', n_epoch);
Te_AUC = zeros(n_epoch,1);
for ei=1:n_epoch
    showprogress(ei, n_epoch);
    nei_lst = nei_info.SubNet_Full{ei};
    eTr = zTr(:, nei_lst); 
    eTe = zTe(:, nei_lst);
    
    switch method_name
        case 'PCA1'
            coeff = pca(eTr, 'NumComponents', 1);
            pred = eTe*coeff(:,1);
        otherwise
            error('Unknown method.');
    end
    
    %% Evaluate prediction
    Te_AUC(ei,1) = getAUC(lTe, pred, 50);
    %[~, ~, ~, auc] = perfcurve(lTe, pred, -1)
end

%% Saving results
sav_name = sprintf('./Result_Files/RES_%s_%s_%s-%s_CVT%02d_Si%02d-Ri%03d.mat', method_name, net_name, tr_name, te_name, cv_ind, Target_Study, Target_Repeat);
fprintf('Saving results in [%s]\n', sav_name);
save(sav_name, 'Te_AUC', 'cv_info', 'nei_info');
end