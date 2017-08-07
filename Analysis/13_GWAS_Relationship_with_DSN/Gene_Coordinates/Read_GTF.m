clc;
clear;

%% Load annotation files
fid = fopen('Gene_Info.tsv', 'r');
gtf_cell = textscan(fid, '%s%s%s%f%f%s%s%s%s', 'HeaderLines', 0, 'Delimiter', '\t', 'CommentStyle', '@', 'ReturnOnError', 0);
if any(~strcmp(gtf_cell{3},'gene')) || ~feof(fid)
	error(); 
end
fclose(fid);

%% Create Gene info dataset
Gene_Crd = [str2double(gtf_cell{1}) gtf_cell{4} gtf_cell{5} strcmp(gtf_cell{7}, '+')*2-1];
Gene_Name = regexp(gtf_cell{9}, '.*gene_name "(.*?)"', 'tokens'); Gene_Name = [Gene_Name{:}]; Gene_Name = [Gene_Name{:}]';
if numel(Gene_Name)~=size(Gene_Crd,1)
	error();
end

%% Filter mitochondria genes
del_ind = isnan(Gene_Crd(:,1));
Gene_Crd(del_ind,:) = [];
Gene_Name(del_ind) = [];

%% Sorting
[~, sid] = sortrows(Gene_Crd(:, 1:2));
Gene_Crd = Gene_Crd(sid,:);
Gene_Name = Gene_Name(sid);

%% Save dataset
save('Gene_Info_vGRCh37.75.mat', 'Gene_Crd', 'Gene_Name');