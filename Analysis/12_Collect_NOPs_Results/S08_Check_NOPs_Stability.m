clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop');
addpath('../../../../Useful_Sample_Codes/FisherExactTest/');

%% Select methods and networks
res_path = './Collected_Results/';
data_lst = {
%     'MRK_CVT01_Lasso_Random-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_Random-P10000_MSN-500.mat'
    
%     'MRK_CVT01_Lasso_I2D-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_I2D-G11748_MSN-500.mat'
    
%     'MRK_CVT01_Lasso_MSigDB-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_MSigDB-G11748_MSN-500.mat'
     
%     'MRK_CVT01_Lasso_HPRD-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_HPRD-G11748_MSN-500.mat'
    
%     'MRK_CVT01_Lasso_KEGG-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_KEGG-G11748_MSN-500.mat'
    
%     'MRK_CVT01_Lasso_STRING-G11748_MSN-500.mat'
    'MRK_CVT01_GLasso_STRING-G11748_MSN-500.mat'
    
%     'MRK_CVT01_Lasso_ACr-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_ACr-P10000_MSN-500.mat'
    
%     'MRK_CVT01_Lasso_AvgSynACr-P10000_MSN-500.mat'
%     'MRK_CVT01_GLasso2_AvgSynACr-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso_AvgSynACr-P10000_MSN-500.mat'
    'MRK_CVT01_GLasso7_AvgSynACr-P10000_MSN-500.mat'
%     'MRK_CVT01_GLasso10_AvgSynACr-P10000_MSN-500.mat'
%     'MRK_CVT01_GLasso20_AvgSynACr-P10000_MSN-500.mat'
        
    'MRK_CVT01_LExAG_None-G11748_MSN-500.mat'
    };
n_data = numel(data_lst);
n_study = 14;
n_rep = 10;

%% Measure stability
X_lbl = {};
Method_FET = {};
for di=1:n_data
    res_name = [res_path data_lst{di}];
    res_info = regexp(data_lst{di}, '_', 'split');
    fprintf('Reading [%s]\n', res_name);
    res_data = load(res_name);
    n_res = size(res_data.Marker_lst,1);
    if any(isnan(res_data.AUC_mat(:))), error(); end
    if n_res~=n_study*n_rep, error(); end
    if ~exist('Gene_Map', 'var')
        Gene_Map = res_data.Gene_Map;
        All_Gene_List = (1:Gene_Map.Count)';
        n_gene = All_Gene_List;
    else
        if ~isequal(Gene_Map, res_data.Gene_Map), error(); end
    end
    
    Top_List = {};
    for si=1:n_study
        is_in = res_data.Marker_SiRi(:,1)==si; % && res_data.Marker_SiRi(:,2)==ri
        if sum(is_in)~=n_rep, error(); end
        Mrk_lst = nonzeros(res_data.Marker_lst(is_in,1:100));
        [Top_Mrk, Top_Freq] = getTop(Mrk_lst, 100);
        
        %if numel(Top_Mrk)<100
        %    Top_Mrk = [Top_Mrk; randi(n_gene, 100-numel(Top_Mrk), 1)];
        %end
        Top_List{di, si} = nonzeros(Top_Mrk);
    end
    
    pval_lst = ones(n_study, n_study);
    for si=1:n_study
        for sj=si+1:n_study
            [~, ~, pval_lst(si,sj)] = getFET(Top_List{di,si}, Top_List{di,sj}, All_Gene_List);
        end
    end
    Method_FET{di,1} = nonzeros(-log10(pval_lst));
    X_lbl{di,1} = sprintf('%s-%s', res_info{3}, res_info{4});
end

%% Plot P-values
close all
figure('Position', [100 100 1500 400]);
hold on
clr_map = jet(n_data); %[AdvancedColormap('cb', n_net/2); AdvancedColormap('mr', n_net/2)];
for di=1:n_data
	box_h = boxplot(Method_FET{di}, 'Position', di, 'Color', clr_map(di,:), 'Widths', 0.5);
    set(box_h, 'LineWidth', 2);
% 	text(di, 250, X_lbl{di}, 'HorizontalAlignment', 'Center', 'VerticalAlignment', 'Top', ...
% 		'FontWeight', 'Bold', 'FontSize', 14);
end
set(gca, 'Xlim', [0 n_data+1], 'Ylim', [0 180], ...
	'XTick', 1:n_data, 'XTicklabel', X_lbl, 'XTicklabelRotation', 10, 'FontWeight', 'Bold');
set(gca,'box','off');
h=findobj(gca,'tag','Outliers');
delete(h);

%% Saving
% output_name = sprintf('./Plots/S08_NOPs_StabilityFET.pdf');
% set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [36 10], 'PaperPosition', [0 0 36 10]);
% print('-dpdf', '-r300', output_name);
