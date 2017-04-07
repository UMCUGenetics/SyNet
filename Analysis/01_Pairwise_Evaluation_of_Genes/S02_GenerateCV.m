function S02_GenerateCV(data_name, cv_type, n_fold, n_rep)
% S02_GenerateCV('SyNet', 1, 10, 20)

GeneExpression_Path = getPath(data_name);
load(GeneExpression_Path, 'Patient_Label', 'Study_Index');
n_sample = numel(Patient_Label);

if exist('Study_Index', 'var')
	n_study = max(Study_Index);
	for ri=1:n_rep
		for si=1:n_study
			cv_obj(si,ri).iTe = Study_Index==si;
			TrnInd = find(Study_Index~=si);
			n_Trn = numel(TrnInd);
			sel_ind = randperm(n_Trn, floor(n_Trn*0.7));
			cv_obj(si,ri).iTr = false(n_sample, 1);
			cv_obj(si,ri).iTr(TrnInd(sel_ind)) = 1;
			if any(cv_obj(si,ri).iTr & cv_obj(si,ri).iTe), error(); end
		end
	end
else
	for ri=1:n_rep
		cv_ind = crossvalind('KFold', Patient_Label, n_fold);
		for fi=1:n_fold
			cv_obj(fi,ri).iTe = cv_ind==fi;
			cv_obj(fi,ri).iTr = cv_ind~=fi;
		end
	end
end

%% Saving
sav_name = sprintf('./CV_Files/CV_%s_CVT%02d.mat', data_name, cv_type);
fprintf('Saving results in [%s]\n', sav_name);
save(sav_name, 'cv_obj', 'Patient_Label');
end
