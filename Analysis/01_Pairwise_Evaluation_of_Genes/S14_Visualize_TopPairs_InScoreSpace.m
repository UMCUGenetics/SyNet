clc;
clear;
close all

%% Initialization
addpath('../../../../Useful_Sample_Codes/getTop/');
addpath('../../../../Useful_Sample_Codes/Advance_Colormap/');
addpath('../../../../Useful_Sample_Codes/OScore/');
addpath('../_Utilities/');
ge_name = 'SyNet';
n_top = 10000;

%% Load Top pairs
net_name = ['./Network_Files/' 'DSN_' ge_name '.mat'];
dsn_info = load(net_name, 'Pair_AUC', 'Gene_Name');
% dsn_info.Gene_Name = dsn_info.Gene_Name(1:50); dsn_info.Pair_AUC = dsn_info.Pair_AUC(1:50,1:50); %###
n_gene = numel(dsn_info.Gene_Name);
n_total = n_gene*(n_gene-1)/2;
Pair_Info = zeros(n_total, 16, 'single');
[Pair_Info(:,1), Pair_Info(:,2)] = find(triu(ones(n_gene), 1));
fprintf('In total [%d] gene pairs exist.\n', n_total);

%% Compute initial axes
Ind_AUC = dsn_info.Pair_AUC(1:n_gene+1:end)';
Pair_Info(:,3) = Ind_AUC(Pair_Info(:,1));
Pair_Info(:,4) = Ind_AUC(Pair_Info(:,2));
Pair_Info(:,5) = max(Pair_Info(:,3:4), [], 2); % Max AUC
pair_ind = sub2ind([n_gene n_gene], Pair_Info(:,1), Pair_Info(:,2));
Pair_Info(:,6) = dsn_info.Pair_AUC(pair_ind); % Combined AUC
Pair_Info(:,7) = Pair_Info(:,6)./Pair_Info(:,5); % Synergy
Pair_Info(:,8) = mean(Pair_Info(:,3:4), 2); % Mean AUC
dsn_info.Pair_AUC = [];
% clear Ind_AUC

%% Add absolute correlation
GE_Path = getPath(ge_name);
fprintf('Loading [%s] expression data\n', GE_Path);
ge_info = load(GE_Path, 'Gene_Expression', 'Gene_Name');
% ge_info.Gene_Name = ge_info.Gene_Name(1:50); ge_info.Gene_Expression = ge_info.Gene_Expression(:,1:50); %###
if ~isequal(ge_info.Gene_Name, dsn_info.Gene_Name), error(); end
Corr_mat = corr(zscore(ge_info.Gene_Expression), 'Type', 'Spearman');
Pair_Info(:,9) = Corr_mat(pair_ind);
Pair_Info(:,10) = abs(Pair_Info(:,9)); % Absolute spearman correlation
clear pair_ind Corr_mat
ge_info.Gene_Expression = [];

%% Normalizing scores
Axes_Name([7 8 10]) = {'Synergy' 'Mean AUC' 'AbsCorr'};
Pair_Info(:,11) = oscore(Pair_Info(:,7)); % Synergy
Pair_Info(:,12) = oscore(Pair_Info(:,8)); % Mean AUC
Pair_Info(:,13) = oscore(Pair_Info(:,10)); % Absolute spearman correlation

%% Normalization step
% Tmp_List = quantilenorm(Pair_Info(:,11:13));

%% Compute final fitness
for ai=11:13
    %Tmp_List(:,ai) = oscore(Tmp_List(:,ai));
    %Tmp_List(isnan(Tmp_List(:,ai)),ai) = 1;
    Pair_Info(:,14) = Pair_Info(:,14) + (1-Pair_Info(:,ai)).^2;
end
clear Tmp_List
Pair_Info(:,14) = -sqrt(Pair_Info(:,14));
Pair_Info(:,15) = oscore(Pair_Info(:,14));
Pair_Info(:,16) = -sqrt(sum((1-Pair_Info(:,11:13)).^2, 2));

%% Sorting
[~, sid] = sort(Pair_Info(:,14), 'Descend');
Pair_Info = Pair_Info(sid, :);

%% Plotting
if 0
    close all
    figure('Visible', 'off');
    ind = [8 7];
    is_top = false(n_total,1); is_top(1:n_top) = 1;
    plot(Pair_Info(~is_top, ind(1)), Pair_Info(~is_top, ind(2)), '.', 'Color', [0.8 0.8 0.8], 'MarkerEdgeColor', [0.8 0.8 0.8], 'MarkerSize', 1);
    hold on
    plot(Pair_Info( is_top, ind(1)), Pair_Info( is_top, ind(2)), '.', 'Color', [0.3 0.8 0.2], 'MarkerEdgeColor', [0.3 0.8 0.2], 'MarkerSize', 4);
    set(gca,'TickDir','out');
    
    [hit_map, center] = hist3([Pair_Info(:,ind(1)) Pair_Info(:,ind(2))], 'Nbins', [50 50]);
    cnt_h = contour(center{1}, center{2}, hit_map', 10, 'ShowText', 'on');
    % cnt_h = contour(center{1}, center{2}, hit_map', 'LineStyle', '--', 'ShowText', 'on', 'LabelSpacing', inf, 'LevelList', logspace(1,7,7));
    clr_map = gray(25); clr_map = clr_map(17:end,:);
    colormap(clr_map);
    
    set(gca, 'FontWeight', 'Bold');
    xlabel(Axes_Name{ind(1)});
    ylabel(Axes_Name{ind(2)});
    title(['Score Space (' ge_name  ')']);
    x_lim = [min(Pair_Info(:,ind(1)))*0.98 max(Pair_Info(:,ind(1)))*1.02];
    y_lim = [min(Pair_Info(:,ind(2)))*0.98 max(Pair_Info(:,ind(2)))*1.02];
    xlim(x_lim);
    ylim(y_lim);
    
    %% Saving the plot
    sav_name = sprintf('./Plots/ScoreSpace_%s_TripleScores_%s-%s.png', ge_name, Axes_Name{ind(1)}, Axes_Name{ind(2)});
    print(gcf, '-dpng', '-r300', sav_name);
end

%% Save the network in tsv file
if 0
    tsv_name = sprintf([net_name(1:end-4) '_TopPairs_N%d.tsv'], n_top);
    fid = fopen(tsv_name, 'w');
    fprintf(fid, 'Source\tTarget\tSynergy\tMean_AUC\tCorr\tAbsCorr\tFitness\tWeight\n');
    for pi=1:n_top
        fprintf(fid, '%s\t%s\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.3f\t%0.6f\n', ge_info.Gene_Name{Pair_Info(pi,1)}, ge_info.Gene_Name{Pair_Info(pi,2)}, Pair_Info(pi,[7 8 9 10 14 15]));
    end
    fclose(fid);
end

%% Save mat file
if 1
    out_mat.PP_Info = Pair_Info(1:n_top*20, :);
    out_mat.NP_Info = Pair_Info(randperm(n_total, n_top*50), :);
    out_mat.Gene_Name = ge_info.Gene_Name;
    out_mat.Ind_AUC = Ind_AUC;
    save(['./Top_Pairs/TopP_' ge_name '.mat'], '-struct', 'out_mat');
    
%     Gene_Name = ge_info.Gene_Name;
%     save(['./Network_Files/AllPairs_' ge_name '.mat'], 'Pair_Info', 'Gene_Name', '-v7.3');
end