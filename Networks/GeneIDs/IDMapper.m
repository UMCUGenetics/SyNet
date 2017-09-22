function [resMap, n_lnk] = IDMapper(mapFilePath, fromCol, toCol, format_str, n_ignore_lines)
if ~exist('format_str', 'var')
    format_str = '%s\t%s\t%s';
end
if ~exist('n_ignore_lines', 'var')
    n_ignore_lines = 0;
end

%% Load file
fid = fopen(mapFilePath, 'r');
res_cell = textscan(fid, format_str, 'MultipleDelimsAsOne', 0, 'Delimiter', '\t', 'HeaderLines', n_ignore_lines, 'ReturnOnError', 0);
fclose(fid);

%% Select columns
fCell = res_cell{:,fromCol};
tCell = res_cell{:,toCol};
n_lnk = numel(fCell);

%% Generating map object
resMap = containers.Map();
for li=1:n_lnk
    if strcmp(fCell{li},'') || strcmp(tCell{li},'')
        continue;
    end
    if resMap.isKey(fCell{li})
        resMap(fCell{li}) = [resMap(fCell{li}) tCell(li)];
    else
        resMap(fCell{li}) = tCell(li);
    end
end
end