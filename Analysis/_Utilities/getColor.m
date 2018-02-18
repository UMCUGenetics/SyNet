function [Method_Color, Preferred_Name] = getColor(Method_Name)
% clr_map = min(hsv(n_res)*1.0, ones(n_res,3));
Preferred_Name = Method_Name;

switch Method_Name
    case {'Random' 'HBGland-SHFL'}
        Method_Color = [0.60 0.60 0.60];
    case 'AvgSynACrNShuff'
        Method_Color = [0.60 0.60 0.60];
        Preferred_Name = 'SyNet-Shuffled';
    case 'HumanInt'
        %Method_Color = [0.63 0.40 0.22]; %old
        Method_Color = [0.62 0.45 0.33];
    case 'BioPlex'
        Method_Color = [0.85 0.60 0.00];
    case 'BioGRID'
        Method_Color = [1.00 0.92 0.00];
    case 'IntAct'
        Method_Color = [0.62 1.00 0.00];
    case 'STRING'
        Method_Color = [0.15 0.80 0.10];
    case 'HBBrain'
        Method_Color = [0.37 0.60 0.70];
        Preferred_Name = 'Brain';
    case 'HBKidney'
        Method_Color = [0.65 0.75 0.14];
        Preferred_Name = 'Kidney';
    case 'HBOvary'
        Method_Color = [0.98 0.48 0.15];
        Preferred_Name = 'Ovary';
    case 'HBGland'
        Method_Color = [1.00 0.20 0.80];
        Preferred_Name = 'Breast';
    case 'HBLympNode'
        %Method_Color = [0.62 0.00 1.00]; % old
        Method_Color = [0.70 0.30 0.96];
        Preferred_Name = 'LymphNode';
    case {'AbsCorr' 'ACr'}
        Method_Color = [0.30 0.50 1.00];
        Preferred_Name = 'Corr';
    case {'SyNet' 'AvgSynACr'}
        Method_Color = [1.00 0.10 0.10];
        Preferred_Name = 'SyNet';
    case {'All genes' 'None'}
        Method_Color = [0.40 0.60 0.40];
        Preferred_Name = 'Baseline';
    case 'STRINGNShuff'
        Method_Color = [0.80 0.80 0.80];
        Preferred_Name = 'STRING-Shuff';
    case {'Syn' 'AvgSyn'}
        Method_Color = [1.00 0.00 0.00];
        case 'HBEpith'
        Method_Color = [0.80 0.40 0.10];
        Preferred_Name = 'Epith';
    case 'HBLiver'
        Method_Color = [0.00 0.90 0.60];
        Preferred_Name = 'Liver';
    case {'All networks'}
        Method_Color = [0.00 0.00 0.00];
    otherwise
        warning('Random color is made for method: [%s]', Method_Name);
        Method_Color = rand(1,3);
end
end