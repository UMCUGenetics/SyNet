addpath(genpath('../../Useful_Sample_Codes/SLEP/'));
addpath(genpath('../../Useful_Sample_Codes/getAUC/'));
addpath(genpath('../../Useful_Sample_Codes/ShowProgress/'));
n_rep = 10*20*1000;

t_ex = 0;
t_mat = 0;
x = rand(3000, 2);
l = (rand(3000, 1)>0.5)*2-1;
z = zscore(x);
for i=1:n_rep
	showprogress(i, n_rep);
	tic;
	B(:,1) = regress(l, z);
	t_mat = t_mat + toc;
	
% 	tic;
% 	B(:,2) = lassoEx(z, l, 'lassoType', 't');
% 	t_ex = t_ex + toc;
end
fprintf('Matlab time = %0.1fs\nEx time = %0.1fs\n', t_mat, t_ex);