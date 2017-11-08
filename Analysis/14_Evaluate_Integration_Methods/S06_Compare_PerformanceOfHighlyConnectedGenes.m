clc;
clear;

%% Initialization
n_study = 14;
n_rep = 5;

%% Infer shuffled genes
net_STRING = load('./NetNei_Files/NetNei_STRING_NN20.mat');
net_Shuffl = load('./NetNei_Files/NetNei_Shf-STRING_NN20.mat');
n_gene = numel(net_STRING.Gene_Name);
rnd_ID = zeros(n_gene, 1);
for gi=1:n_gene
    rnd_ID(gi) = net_Shuffl.SubNet_Full{gi}(1);
end

%% Get gene degree
for gi=1:n_gene
    shf_ind = net_STRING.SubNet_Full{gi};
    str_ind = net_Shuffl.SubNet_Full{gi};
    
    if ~isequal(str_ind, rnd_ID(shf_ind))
        error();
    end
end
str_degree = cellfun('length', net_STRING.SubNet_Full);
shf_degree = cellfun('length', net_Shuffl.SubNet_Full);
if ~isequal(str_degree, shf_degree)
    exist();
end

%% Collection of data
ind_auc = zeros(n_study, n_rep, n_gene);
for si=1:n_study
    for ri=1:n_rep
        res_name = sprintf('ResIND_SyNet-SyNet_CVT50_Si%02d-Ri%03d.mat', si, ri);
        fprintf('Result: %s\n', res_name);
        
        res_data = load(['./Result_Files/' res_name]);
        ind_auc(si, ri, :) = res_data.Gene_TrAUC;
    end
end
ind_auc = squeeze(mean(mean(ind_auc, 1),2));

% %% Sort genes based on degree
% [str_val, str_sid] = sort(str_degree, 'Descend');
% [shf_val, shf_sid] = sort(shf_degree, 'Descend');
% 
% %% Measure performance of genes
% str_auc = ind_auc(str_sid);
% shf_auc = ind_auc(shf_sid);

%% Plotting
close all
figure('Position', [100 100 1500 700]);
hold on
size_set = floor(linspace(20, 1, 7))';
n_step = size(size_set, 1);
size_class = [size_set(1:end-1) size_set(2:end)-1];
clr_map = [
    1 0 0;
    0 0 1;
    ];
for si=1:n_step
    has_ol = str_degree>=size_class(si,1) & str_degree<=size_class(si,2);
    boxplot(ind_auc(has_ol), 'Positions', si-0.2, 'Color', clr_map(1,:));
    
    has_ol = shf_degree>=size_class(si,1) & shf_degree<=size_class(si,2);
    boxplot(ind_auc(has_ol), 'Positions', si+0.2, 'Color', clr_map(2,:));
end
xlim([0 n_step+1]);

%% Saving the plot
output_name = sprintf('./Plots/S06_PerformanceComparison_STRING_PerDegree.pdf');
set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [16 4], 'PaperPosition', [0 0 16 4]);
print('-dpdf', '-r300', output_name);

