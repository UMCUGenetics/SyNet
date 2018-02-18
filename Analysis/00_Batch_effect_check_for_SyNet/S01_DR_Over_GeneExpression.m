clc;
clear;

%% Initialization
data_lst = {
    '../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat'
    '../../Gene_Expression_Datasets/SyNet/SyNet_Normalized_Par.mat'
    };
data_name = {'BatchCorrected' 'Normalized'};
n_data = numel(data_lst);
mrk_lst = {'o', '+', '*', 'v', '.', 'x', 's', 'd', '^', '>', '<', 'p', 'h'};
n_mrk = numel(mrk_lst);

%% Main loop
for di=1:n_data
    
    %% Load gene expression
    ge_path = data_lst{di};
    fprintf('Loading gene expression in [%s]\n', ge_path);
    load(ge_path, 'Gene_Expression', 'Patient_Info');
    
    if di==1
        invalid_survival = isnan(Patient_Info.Prognostic_Status);
        Patient_Info(invalid_survival, :) = [];
        Gene_Expression(invalid_survival, :) = [];
        is_in = ismember(Patient_Info.Source_Study, {'METABRIC' 'TCGA'});
        is_in(1:1616) = 1;
        Patient_Info = Patient_Info(is_in,:);
        Gene_Expression = Gene_Expression(is_in,:);
    end
    
    %[Study_Name, ~, Study_Index] = unique(strcat(Patient_Info.Source_Study, ';', Patient_Info.StudyName), 'Stable');
    [Study_Name, ~, Study_Index] = unique(Patient_Info.StudyName, 'Stable');
    Patient_Label = Patient_Info.Prognostic_Status;
    zData = zscore(Gene_Expression);
    if any(isnan(zData(:))), error(); end
    n_study = max(Study_Index);
    
    %% Plot class distribution
    if 0
        cls_dist = zeros(n_study, 2);
        for si=1:n_study
            cls_dist(si,1) = sum(Study_Index==si & Patient_Label==0);
            cls_dist(si,2) = sum(Study_Index==si & Patient_Label==1);
        end
        close all
        bar(cls_dist, 'Grouped');
        colormap([0.3 0.3 1; 1 0.3 0.3]);
        ylabel('# Patients');
        set(gca, 'FontWeight', 'Bold', 'XTick', 1:n_study, 'XTickLabel', Study_Name, 'XTickLabelRotation', 25, 'XLim', [0 n_study+1]);
        legend({'Good' 'Poor'}, 'fontsize', 12, 'fontweight', 'bold', 'Location', 'NorthWest');
        set(gcf, 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [60 10], 'PaperPosition', [0 0 60 10], 'PaperUnits', 'Centimeter');
        sav_name = sprintf('./Plots/ClassFrequency_%s.pdf', data_name{di});
        print('-dpdf', '-r300', sav_name);
    end
    
    %% Reduce dimention
    perp_lst = [20]; % 5 10 20 40 70 110 200 400
    for pi=1:numel(perp_lst)
        perplexity = perp_lst(pi);
        fprintf('Perplexity set to: %d\n', perplexity);
        initial_dims = 30;
        no_dims = 2;
        pData = fast_tsne(zData, no_dims, initial_dims, perplexity, 0.5, 500, 1);
        
        %% Plot Outcome
        if 0
            close all
            clr_map = [0 0 1; 1 0 0];
            figure('Visible', 'off');
            hold on
            for ci=1:2
                plot(pData(Patient_Label==ci-1,1), pData(Patient_Label==ci-1,2), mrk_lst{ci}, 'markeredgecolor', clr_map(ci,:), 'markersize', 3, 'LineWidth', 0.4);
            end
            legend({'Good' 'Poor'}, 'fontsize', 12, 'fontweight', 'bold');
            sav_name = sprintf('./Plots/GeneExpression_%s_PatientLabel_Perp%03d.pdf', data_name{di}, perplexity);
            print(gcf, '-dpdf', '-r300', sav_name);
        end
        
        %% Plot Studies
        if 1
            close all
            clr_map = jet(n_study) * 0.8;
            %clr_map = lines(n_study) * 0.8;
            clr_map = [winter(n_study-2); [97, 26, 132]/255; 1 0 0];
            figure('Position', [50 50 1000 700], 'Visible', 'on');
            hold on
            mrk_h = [];
            for si=1:n_study
                plot(pData(Study_Index==si,1), pData(Study_Index==si,2), mrk_lst{mod(si-1,n_mrk)+1}, 'markerfacecolor', 'none', 'markeredgecolor', clr_map(si,:), 'markersize', 3);
                mrk_h(si) = plot(nan, nan, mrk_lst{mod(si-1,n_mrk)+1}, 'markerfacecolor', 'none', 'markeredgecolor', clr_map(si,:), 'markersize', 7, 'LineWidth', 1.5);
            end
            legend(mrk_h, Study_Name, 'fontsize', 8, 'fontweight', 'bold', 'Location', 'EastOutside');
            set(gcf, 'PaperUnits', 'Inches', 'PaperOrientation', 'landscape', 'PaperPositionMode','auto', 'PaperSize', [12 10], 'PaperPosition', [0 0 12 10]);
            sav_name = sprintf('./Plots/GeneExpression_%s_Study_Perp%03d.pdf', data_name{di}, perplexity);
            print(gcf, '-dpdf', '-r300', sav_name);
        end
    end
end


