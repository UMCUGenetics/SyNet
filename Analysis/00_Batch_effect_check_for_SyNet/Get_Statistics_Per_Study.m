clc;
clear

%% Initialization
data_lst = {
    '../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat'
    '../../Gene_Expression_Datasets/HaibeKains/HaibeKains_Combined.mat'
    };
n_data = numel(data_lst);

%% Main loop
Study_Info = {};
step = 1;
for di=1:n_data
    fprintf('Reading [%s] ...\n', data_lst{di});
    load(data_lst{di}, 'Patient_Info');
    Study_List = Patient_Info.StudyName;
    if di==1
        % Surv_Time = Patient_Info.SurvivalTime;
        Prognostic_Status = Patient_Info.Prognostic_Status;
    else
        Surv_Time = str2double([Patient_Info.DMFSTime Patient_Info.RFSTime Patient_Info.OSTime]);
        for ci=2:size(Surv_Time,2)
            is_nan = isnan(Surv_Time(:,1));
            if ~any(is_nan), break; end
            Surv_Time(is_nan,1) = Surv_Time(is_nan,ci);
        end
        Prognostic_Status = double(Surv_Time(:,1)<=1825);
        Prognostic_Status(isnan(Surv_Time(:,1))) = nan;
    end
    
    %if any(~ismember(Prognostic_Status, [0 1 nan])), error(); end
    for si_name = unique(Study_List, 'Stable')'
        fprintf('\t*** Calculating for [%s]\n', si_name{1});
        is_in = strcmp(Study_List, si_name{1});
        Study_Info{step, 1} = si_name{1}; % Study name
        Study_Info{step, 2} = sum(is_in); % # Samples
        Study_Info{step, 3} = sum(Prognostic_Status(is_in)==1); % # Poor outcome
        Study_Info{step, 4} = sum(Prognostic_Status(is_in)==0); % # Good outcome
        Study_Info{step, 5} = sum(~isnan(Prognostic_Status(is_in))); % # Useable sample
        Study_Info{step, 6} = sum( isnan(Prognostic_Status(is_in))); % # NaN survival
        step = step + 1;
    end
end

