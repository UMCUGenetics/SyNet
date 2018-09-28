
fn_lst = dir('./Results_Files/*SyNet*HubGL5*.mat');
n_file = numel(fn_lst);
n_gene = zeros(n_file, 1);
for fi=1:n_file
	fname = ['./Results_Files/' fn_lst(fi).name];
	fprintf('Loading %s: ', fname);
	data = load(fname);
	n_gene(fi) = numel(data.Gene_Name);
	fprintf('%d\n', n_gene(fi));
end
n_gene
