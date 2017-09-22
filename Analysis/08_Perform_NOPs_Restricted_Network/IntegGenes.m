function [cTr, cTe] = IntegGenes(Nei_grp, zTr, lTr, zTe, Study_Index, Int_Type)

%% Prepare data
% if ~exist('Int_Type', 'var'), Int_Type='Mean'; end
% fprintf('Integrating genes using [%s] method ...\n', Int_Type);
n_zTr = size(zTr, 1);
n_zTe = size(zTe, 1);
n_grp = numel(Nei_grp);
[~, ~, Fold_Index] = unique(Study_Index, 'Stable');
n_Fold = max(Fold_Index);

%% Combination
cTr = zeros(n_zTr, n_grp);
cTe = zeros(n_zTe, n_grp);
warning('off');
for gi=1:n_grp
	gene_set = Nei_grp{gi};
	set_size = numel(gene_set);
	
	if strcmpi(Int_Type, 'RI')
		B_fold = zeros(set_size, n_Fold);
		for fi=1:n_Fold
			B_fold(:, fi) = regress(lTr(Fold_Index==fi), zTr(Fold_Index==fi, gene_set));
		end
		B_opt = median(B_fold, 2);
		cTr(:, gi) = zTr(:, gene_set) * B_opt;
		cTe(:, gi) = zTe(:, gene_set) * B_opt;
	else
		cTr(:, gi) = mean(zTr(:, gene_set), 2);
		cTe(:, gi) = mean(zTe(:, gene_set), 2);
	end	
end
warning('on');
end