function S02_GenerateCV(data_name, cv_type, n_rep)
% S02_GenerateCV('SyNet', 2, 5)
%% Initialization
addpath('../_Utilities/');

GeneExpression_Path = getPath(data_name);
load(GeneExpression_Path, 'Patient_Label', 'Study_Index');
n_sample = numel(Patient_Label);
n_study = max(Study_Index);

for vi=1:n_study
	step = 1;
	for si=1:n_study
		if vi==si, continue; end
		for ri=1:n_rep
			cv_obj(vi,step).iTv = Study_Index==vi;
			cv_obj(vi,step).iTe = Study_Index~=vi & Study_Index==si;
			TrnInd =         find(Study_Index~=vi & Study_Index~=si);
			n_Trn = numel(TrnInd);
			cv_obj(vi,step).iTr = false(n_sample, 1);
			sel_ind = randperm(n_Trn, floor(n_Trn*0.7));
			cv_obj(vi,step).iTr(TrnInd(sel_ind)) = 1;
			if sum([cv_obj(vi,step).iTr cv_obj(vi,step).iTe cv_obj(vi,step).iTv],2)>1, error(); end
			step = step + 1;
		end
	end
end

%% Saving
sav_name = sprintf('./CV_Files/CV_%s_CVT%02d.mat', data_name, cv_type);
fprintf('Saving results in [%s]\n', sav_name);
save(sav_name, 'cv_obj', 'Patient_Label', 'Study_Index');
end
