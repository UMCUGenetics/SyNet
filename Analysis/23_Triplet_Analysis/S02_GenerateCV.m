function S02_GenerateCV()

%% Initialization
addpath('../_Utilities/');
n_rep = 5;
cv_type = 'AcrossStudies';
data_name = 'SyNet';

%% Loading data indices
GeneExpression_Path = getPath(data_name);
load(GeneExpression_Path, 'Patient_Label', 'Study_Index');
n_sample = numel(Patient_Label);
n_study = max(Study_Index);

%% Main loop
for si=1:n_study
    for ri=1:n_rep
        cv_obj(si,ri).iTe = Study_Index==si;
        TrnInd = find(Study_Index~=si);
        n_Trn = numel(TrnInd);
        cv_obj(si,ri).iTr = false(n_sample, 1);
        sel_ind = randperm(n_Trn, floor(n_Trn*0.7));
        cv_obj(si,ri).iTr(TrnInd(sel_ind)) = 1;
        if sum([cv_obj(si,ri).iTr cv_obj(si,ri).iTe],2)>1, error(); end
    end
end

%% Saving
sav_name = sprintf('./CV_Files/CV_%s_CVT-%s.mat', data_name, cv_type);
fprintf('Saving results in [%s]\n', sav_name);
save(sav_name, 'cv_obj', 'Patient_Label', 'Study_Index');
end
