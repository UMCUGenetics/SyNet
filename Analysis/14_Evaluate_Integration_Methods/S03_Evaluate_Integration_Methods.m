function S03_Evaluate_Integration_Methods(Target_Study, Target_Repeat, method_name, net_name, MIN_SUBNET_SIZE)
%% Run
%{
for mi in Avg Std DA2 Reg PCA1 DPCA DA2NoRem Rnd; do
for ri in `seq 1 5`; do
for si in `seq 1 14`; do
PARAM="$si,$ri,'$mi','STRING_NN20',20"; 
sbatch --job-name=IE-$PARAM --output=Logs/NE-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=5GB --time=04:00:00 --ntasks=1 --cpus-per-task=1 run_Matlab.sh S03_Evaluate_Integration_Methods "$PARAM";
done;
done;
read -p "Press a key" -t 1800;
done
%}


%% Initialization
clc
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
addpath('../../../../Useful_Sample_Codes/getAUC/');
addpath('../_Utilities/');
addpath(genpath('../../../../Useful_Sample_Codes/SLEP'));
tr_name = 'SyNet';
te_name = 'SyNet';
cv_ind = 50;
if ismac
    method_name = 'Rnd';
    net_name = 'STRING_NN05';
    Target_Study = 3;
    Target_Repeat = 2;
    MIN_SUBNET_SIZE = 5;
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
n_gene = size(zTr, 2);

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
nei_info.nei_name = nei_name;

%% Filter subnetworks
Grp_Size = cellfun('length', nei_info.SubNet_Full);
nei_info.SubNet_Full(Grp_Size<MIN_SUBNET_SIZE+1) = [];
n_epoch = numel(nei_info.SubNet_Full);

%% Get individual AUC
resind_name = sprintf('./Result_Files/ResIND_%s-%s_CVT%02d_Si%02d-Ri%03d.mat', tr_name, te_name, cv_ind, Target_Study, Target_Repeat);
if exist(resind_name, 'file')
    fprintf('Individual AUCs are found in [%s]\n', resind_name);
    load(resind_name, 'Gene_TrAUC');
else
    fprintf('Evaluating AUC of [%d] individual genes:\n', n_gene);
    Gene_TrAUC = zeros(n_gene,1);
    for gi=1:n_gene
        showprogress(gi, n_gene);
        Gene_TrAUC(gi,1) = getAUC(lTr, zTr(:,gi), 50);
    end
    save(resind_name, 'Gene_TrAUC');
end

%% Main loop
fprintf('Comparing performance for [%d] neighbor sets:\n', n_epoch);
Te_AUC = zeros(n_epoch,1);
Ind_AUC = zeros(n_epoch,1);
Used_IND = cell(n_epoch,1);
for ei=1:n_epoch
    showprogress(ei, n_epoch);
    nei_lst = nei_info.SubNet_Full{ei};
    eTr = zTr(:, nei_lst);
    eTe = zTe(:, nei_lst);
    grp_size = numel(nei_lst);
    Used_IND{ei} = nei_lst;
    
    switch method_name
        case 'Avg'
            pred = mean(eTe, 2);
        case 'Std'
            pred = std(eTe, 0, 2);
        case 'DA2'
            gene_dir = corr(eTr, lTr, 'Type', 'Spearman');
            pTe = eTe;
            for pi=1:grp_size
                if gene_dir(pi)<0
                    pTe(:, pi) = -eTe(:, pi);
                end
            end
            val_ind = abs(gene_dir)>0.02;
            if any(val_ind)
                pred = mean(pTe(:, val_ind), 2);
            else
                pred = mean(pTe, 2);
            end
        case 'DA2NoRem'
            gene_dir = corr(eTr, lTr, 'Type', 'Spearman');
            pTe = eTe;
            for pi=1:grp_size
                if gene_dir(pi)<0
                    pTe(:, pi) = -eTe(:, pi);
                end
            end
            pred = mean(pTe, 2);   
        case 'Reg'
            B = regress(lTr, eTr);
            pred = eTe * B;
        case 'PCA1'
            coeff = pca(eTr, 'NumComponents', 1);
            pred = eTe*coeff(:,1);
        case 'DPCA'
            coeff = pca(eTr);
            pTr = eTr*coeff;
            pc_auc = zeros(1, grp_size);
            for pi=1:grp_size
                pc_auc(pi) = getAUC(lTr, pTr(:,pi), 50);
            end
            [~, tp_ind] = max(pc_auc);
            pred = eTe*coeff(:, tp_ind);
        case 'Rnd'
            rnd_ind = randi(grp_size, 1);
            pred = eTe(:, rnd_ind);
        otherwise
            error('Unknown method.');
    end
    
    %% Evaluate prediction
    Te_AUC(ei,1) = getAUC(lTe, pred, 50);
    %[~, ~, ~, auc] = perfcurve(lTe, pred, -1)
    [~, tg_ind] = max(Gene_TrAUC(nei_lst));
    Ind_AUC(ei) = getAUC(lTe, eTe(:,tg_ind), 50);
end

%% Saving results
sav_name = sprintf('./Result_Files/RES_%s_%s_%s-%s_CVT%02d_Si%02d-Ri%03d.mat', method_name, net_name, tr_name, te_name, cv_ind, Target_Study, Target_Repeat);
fprintf('Saving results in [%s]\n', sav_name);
save(sav_name, 'Te_AUC', 'cv_info', 'nei_info', 'Ind_AUC', 'Gene_TrAUC', 'Used_IND');
end