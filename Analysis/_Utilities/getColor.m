function [Method_Color, Preferred_Name] = getColor(Method_Name)
% clr_map = min(hsv(n_res)*1.0, ones(n_res,3));
Preferred_Name = Method_Name;

switch Method_Name
    case 'Random'
        Method_Color = [1.00 0.00 0.00];
        Preferred_Name = 'STRING-Suff';
    case 'BioPlex'
        Method_Color = [1.00 0.50 0.00];
    case 'BioGRID'
        Method_Color = [1.00 1.00 0.00];
    case 'IntAct'
        Method_Color = [0.50 1.00 0.00];
    case 'STRING'
        Method_Color = [0.00 1.00 0.00];
    case 'HBGland'
        Method_Color = [0.00 1.00 0.50];
        Preferred_Name = 'Gland';
    case 'HBLympNode'
        Method_Color = [0.00 1.00 1.00];
        Preferred_Name = 'LympNode';
    case 'HBEpith'
        Method_Color = [0.00 0.50 1.00];
        Preferred_Name = 'Epith';
    case {'AbsCorr' 'ACr'}
        Method_Color = [0.00 0.00 1.00];
        Preferred_Name = 'AbsCorr';
    case {'tTest'}
        Method_Color = [0.50 0.00 1.00];
    case {'All genes' 'None'}
        Method_Color = [1.00 0.00 1.00];
        Preferred_Name = 'All genes';
    case {'SyNet' 'AvgSynACr'}
        Method_Color = [1.00 0.00 0.50];
        Preferred_Name = 'SyNet';
    case {'Syn' 'AvgSyn'}
        Method_Color = [1.00 0.50 0.75];
    otherwise
        error('Unknown method: [%s]', Method_Name);
end
end