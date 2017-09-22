
%% Initialization
addpath('../01_Pairwise_Evaluation_of_Genes/');
cv_path = './CV_Files/';
cv_type = 1;
geTr_name = 'SyNet';
geTe_name = 'SyNet';

%% Load Training and Test sets
geTr_path = getPath(geTr_name);
fprintf('Training is set to [%s]\n', geTr_path);
tr_info = load(geTr_path, 'Gene_Name', 'Patient_Label');

geTe_path = getPath(geTe_name);
fprintf('Testing is set to [%s]\n', geTe_path);
te_info = load(geTe_path, 'Gene_Name', 'Patient_Label', 'Study_Name', 'Study_Index');

%% Generate CVs
cv_info = load('../01_Pairwise_Evaluation_of_Genes/CV_Files/CV_SyNet_CVT01.mat');
if ~isequal(te_info.Patient_Label, cv_info.Patient_Label) || ~isequal(tr_info.Patient_Label, cv_info.Patient_Label), error(); end
[n_study, n_rep] = size(cv_info.cv_obj);
for si=1:n_study
	for ri=1:n_rep
		cv_item = cv_info.cv_obj(si, ri);
		% rng shuffle % creates a different seed each time, but its not nessasary as seed is set in the job submission
		cv_id = [datestr(now, 'yymmdd') sprintf('Si%02dRi%02d', si, ri)];
		cv_name = sprintf([cv_path 'CID-%s_CVT%02d_ETr-%s_ETe-%s.mat'], cv_id, cv_type, geTr_name, geTe_name);
		if exist(cv_name, 'file')
			error('CV set [%s] is already generated.', cv_name);
		end
		fprintf('Generating CV set [%s].\n', cv_name);
		
		%% Generate Train
		tr_info.CVInd = cv_item.iTr;
		tr_info.GEPath = geTr_path;
		tr_info.GEName = geTr_name;
		tr_info.iCvPar = cvpartition(tr_info.Patient_Label(tr_info.CVInd), 'Kfold', 5);
		fprintf('From total samples of [%d], [%d] are selected.\n', numel(tr_info.CVInd), sum(tr_info.CVInd));
		
		%% Generate Test
		te_info.CVInd = cv_item.iTe;
		te_info.GEPath = geTe_path;
		te_info.GEName = geTe_name;
		te_info.SrcStudy = unique(te_info.Study_Name(te_info.Study_Index(cv_item.iTe)));
		fprintf('From total samples of [%d], [%d] are selected.\n', numel(te_info.CVInd), sum(te_info.CVInd));
		
		%% Save CV info
		if isequal(geTr_name, geTe_name) && any(tr_info.CVInd & te_info.CVInd)
			error();
		end
		save(cv_name, 'tr_info', 'te_info');
		fprintf('CV set is saved in [%s].\n', cv_name);
	end
end
