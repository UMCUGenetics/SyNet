function S03_PerformPWR(batch_be, step_size)
%{
for bi in `seq 1 20 320`; do
PARAM="$bi,20"; 
echo sbatch --exclude=maxwell --job-name=TRC-$PARAM --output=Logs/TRC-$PARAM.%J_%a-%N.out --partition=general --qos=short --mem=5GB --time=04:00:00 --ntasks=1 --cpus=1 --cpus-per-task=1 run_Matlab.sh S03_PerformPWR "$PARAM";
done
%}
if ismac
    batch_be = 10;
    batch_en = 13;
end

%% Initialization
addpath('../../../../Useful_Sample_Codes/getAUC/');
addpath('../../../../Useful_Sample_Codes/fastAUC/');
addpath('../../../../Useful_Sample_Codes/ShowProgress/');
addpath('../_Utilities/');
ge_name = 'SyNet';
trc_path = ['./TRC_Files/' ge_name '/'];
cv_path = './CV_Files/';
Cmb_List = {1 2 3 [1,2] [1,3] [2,3] [1,2,3]};
n_cmb = numel(Cmb_List);
n_candid = 314;
batch_en = batch_be + step_size;
fprintf('/// Running triplet evaluation for [%s], we will test [%d-%d] pairs.\n', ge_name, batch_be, batch_en);

%% Load triplets according to reference gene set
if batch_en > n_candid, batch_en = n_candid; end
load('./Gene_List/Reference_GList.mat', 'Ref_GeneIndex', 'Ref_GeneName');

%% Candid selection
% if ~isequal(data.Gene_Name(Ref_GeneIndex), Ref_GeneName), error(); end
pair_lst = combnk(batch_be:batch_en, 2);
Triplets_Index = [repmat(pair_lst, n_candid, 1) repelem(1:n_candid, size(pair_lst,1))'];
n_triplet = size(Triplets_Index, 1);
fprintf('In total [%d] gene triplets exist.\n', n_triplet);
clear pair_lst

%% Load data
GeneExpression_Path = getPath(ge_name);
fprintf('Loading [%s]\n', GeneExpression_Path);
data = load(GeneExpression_Path, 'Gene_Expression', 'Patient_Label', 'Gene_Name');
zData = zscore(data.Gene_Expression);
Patient_Label = double(data.Patient_Label);
[n_sample, n_gene] = size(zData);
fprintf('Data loaded [n_sample=%d, n_gene=%d]...\n', n_sample, n_gene);
Data_Info.Gene_Name = data.Gene_Name;
Data_Info.Patient_Label = data.Patient_Label;
clear data

%% Load CV info
cv_name = [cv_path 'CV_' ge_name '_CVT-AcrossStudies.mat'];
CV_Info = load(cv_name);
if ~isequal(CV_Info.Patient_Label, Patient_Label), error(); end
cv_obj = CV_Info.cv_obj;
[n_fold, n_rep] = size(cv_obj);
fprintf('Loading CV info from [%s], [%d] folds and [%d] repeats are loaded.\n', cv_name, n_fold, n_rep);

%% Main loop
fprintf('Triplet comparison started at: %s\n', datetime);
fprintf('Evaluating triplets:\n');
Triplet_AUC = zeros(n_triplet, 3+n_cmb);
warning off
for ti=1:n_triplet
    % showprogress(ti, n_selt, 10, '%0.0f%%\n');
    if mod(ti, 1000)==0
        fprintf('Currently at [%06d/%06d] for genes [%04d/%04d/%04d]\n', ti, n_triplet, Triplets_Index(ti,:));
    end
    gene_set = Ref_GeneIndex(Triplets_Index(ti,:));
    auc_cmb = zeros(1, n_cmb);
    for ci=1:n_cmb
        auc_mat = zeros(n_fold, n_rep);
        for fi=1:n_fold
            for ri=1:n_rep
                %fprintf('Rep [%d], fold [%d]\n', ri, fi);
                iTr = cv_obj(fi, ri).iTr;
                iTe = cv_obj(fi, ri).iTe;
                zTr = zData(iTr, gene_set(Cmb_List{ci}));
                lTr = Patient_Label(iTr);
                zTe = zData(iTe, gene_set(Cmb_List{ci}));
                lTe = Patient_Label(iTe);
                
                B = regress(lTr, zTr);
                pred = zTe * B;
                
                % auc = getAUC(lTe, pred);
                % [~,~,~,auc] = perfcurve(lTe, pred, 1)
                tmp = fastAUC(lTe, pred, 1);
                auc_mat(fi, ri) = max([1-tmp tmp]);
            end
        end
        auc_cmb(ci) = mean(mean(auc_mat));
    end
    Triplet_AUC(ti, :) = [Ref_GeneIndex(Triplets_Index(ti,:))' auc_cmb];
end
fprintf('Triplet comparison finished at: %s\n', datetime);
warning on

%% Compute triple synergy
Triplet_AUC(:,11) = Triplet_AUC(:,10) ./ max(Triplet_AUC(:, 4:9),[],2);
[~, sind] = sort(Triplet_AUC(:,11), 'Descend');
Triplet_AUC = Triplet_AUC(sind, :);

%% Saving
sav_name = sprintf('%sTC_%s_%06d-%06d.mat', trc_path, ge_name, batch_be, batch_en);
fprintf('Saving result in [%s]', sav_name);
save(sav_name, 'Triplet_AUC', 'cv_name', 'Cmb_List', 'CV_Info', 'Data_Info', 'Ref_GeneIndex', 'Ref_GeneName');
fprintf('Process finished at: %s\n', datetime);
end

