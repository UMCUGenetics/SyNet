function pv_mat = ttest2Ex(X, Label, varargin)
n_feat = size(X, 2);
Class = unique(Label, 'Stable');
if numel(Class)~=2, error('Label does not contain 2-classes'); end

%% Calculating p-values
pv_mat = zeros(n_feat, 1);
for fi=1:n_feat
    in = Label==Class(1);
    [~, pv_mat(fi)] = ttest2(X(in,fi), X(~in,fi), varargin{:});
end
end