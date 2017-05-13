% wget http://string-db.org/download/protein.links.detailed.v10/9606.protein.links.detailed.v10.txt.gz

% [No Need Anymore]
% zcat 9606.protein.links.detailed.v10.txt.gz | awk '($6!="0"){print($0)}' | sort -k6,6 > STRING_CoExpressionOnly.txt
% after that, bring the header to top

%% Initialization
clc;
clear;
addpath('../../../../Useful_Sample_Codes/GeneID_Coversion/');

%% Load file
fid = fopen('./9606.protein.links.detailed.v10.txt', 'r'); %
Header_List = regexp(fgetl(fid), ' ', 'Split')';
Lnk_lst = textscan(fid, '%s%s%d%d%d%d%d%d%d%d', 'HeaderLines', 0, 'Delimiter', ' ', 'CommentStyle', '@', 'ReturnOnError', 0);
if ~feof(fid)
	error();
end

%% Load the mapper
addpath('../GeneIDs/');
pMap = IDMapper('../GeneIDs/ENSP-ENST-ENSG-HGNC.txt', 1, 4, '%s\t%s\t%s\t%s', 1);

%% ID conversion
fprintf('Converting IDs.\n');
n_line = numel(Lnk_lst{1});
for li=1:n_line
	if pMap.isKey(Lnk_lst{1}{li}(6:end))
		gSrc = pMap(Lnk_lst{1}{li}(6:end));
		Lnk_lst{1}{li} = gSrc{1};
	end
	if pMap.isKey(Lnk_lst{2}{li}(6:end))
		gTar = pMap(Lnk_lst{2}{li}(6:end));
		Lnk_lst{2}{li} = gTar{1};
	end
end

%% Main loop
for hi=3:numel(Header_List)
	fprintf('Generating links for [%s] ... \n', Header_List{hi});
	src_set = Lnk_lst{1};
	tar_set = Lnk_lst{2};
	scr_set = Lnk_lst{hi};
	
	%% Sorting
	[val, sid] = sort(scr_set, 'Descend');
	sid(val==0) = [];
	src_set = src_set(sid);
	tar_set = tar_set(sid);
	scr_set = scr_set(sid);
	n_line = numel(scr_set);
	
	%% Saving
	fid = fopen(sprintf('STRING_%s_GN_WithWeights.txt', Header_List{hi}), 'w');
	for li=1:n_line
		fprintf(fid, '%s\t%s\t%d\n', src_set{li}, tar_set{li}, scr_set(li));
	end
	fclose(fid);
end
fprintf('Process finished.\n');
