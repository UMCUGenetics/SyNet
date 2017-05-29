clc;
clear;

%% Initialization
result_path = '../08_Perform_NOPs_Restricted_Network/Results_Files/';
cv_ind = 1;
method_lst = {'TNMC' 'TLEx' 'Lasso' 'CFGLasso'};
n_met = numel(method_lst);
net_lst = {
	'AvgSynCrr-T10000','Corr-T10000','STRING-T10000','KEGG-T10000','Random-T10000'
	};
n_net = numel(net_lst);
n_study = 14;
n_rep = 10;
MAX_N_SN = 500;
n_TopMarker = 500;
sav_path = './Collected_Markers/';
[~,~] = mkdir(sav_path);

%% Generate gene map
GE_Data = load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name', 'Study_Name');
n_gene = numel(GE_Data.Gene_Name);
global Gene_Map
Gene_Map = containers.Map(GE_Data.Gene_Name, 1:n_gene);

%% Main loop
Result_Marker = zeros(0, 5 + n_TopMarker);
% 1:Method_Index, 2:Network_Index, 3:Study_Index, 4:Repeat_Index, 5:Test_AUC
Result_AUC = nan(n_met, n_net, n_study, n_rep);
for mi=1:n_met
	for ni=1:n_net
		sav_name = sprintf([sav_path 'MRK_CVT%02d_%s_%s_MSN-%03d.mat'], cv_ind, method_lst{mi}, net_lst{ni}, MAX_N_SN);
		if exist(sav_name, 'file')
			fprintf('/// Loading from [%s]\n', sav_name);
			res_data = load(sav_name);
			if ~isequal(res_data.Gene_Map, Gene_Map), error(); end
			Part_AUC = res_data.Part_AUC;
			Part_Marker = res_data.Part_Marker;
		else
			Part_AUC = zeros(n_study, n_rep);
			Part_Marker = [];
			step = 1;
			for si=1:n_study
				for ri=1:n_rep
					res_name = sprintf('%sDID_CVT%02d_Si%02d-Ri%03d_%s_*_MSN-%03d_MTN-%s.mat', result_path, cv_ind, si, ri, net_lst{ni}, MAX_N_SN, method_lst{mi});
					file_info = dir(res_name);
					if numel(file_info)==0
						fprintf('************ Ignoring [%60s]\n', res_name);
						continue;
					elseif numel(file_info)>1
						fprintf('+++ Warning, multiple files found for [%60s]\n', res_name);
						[~, sind] = sort([file_info(:).datenum], 'Descend');
						file_info = file_info(sind);
						for i=2:numel(file_info)
							delete([result_path file_info(i).name]);
						end
					end
					fprintf('Reading from: [%s]\n', file_info(1).name);
					res_data = load([result_path file_info(1).name]);
					switch method_lst{mi}
						case 'TNMC'
							SubNet_Score = res_data.SubNet_Score;
						case 'TLEx'
							SubNet_Score = abs(res_data.B(:,res_data.fit.IndexMinMSE));
						case 'Lasso'
							SubNet_Score = abs(res_data.B(:,res_data.fit.IndexMinMSE));
							res_data.SubNet_List = num2cell(1:numel(res_data.Gene_Name))';
						case 'CFGLasso'
							Group_Index = res_data.fit.Options.ind; %res_data.SubNet_GInd;
							n_snet = size(Group_Index,2);
							SubNet_Score = zeros(n_snet, 1);
							for gi=1:n_snet
								set_ind = Group_Index(1,gi):Group_Index(2,gi);
								SubNet_Score(gi) = max(abs(res_data.B(set_ind, res_data.fit.IndexMinMSE)));
							end
						case 'LExAG'
							error('Check code');
							SubNet_Score = abs(res_data.B(:,res_data.fit.IndexMinMSE));
							res_data.SubNet_List = num2cell(1:n_gene)';
							res_data.Gene_Name = GE_Data.Gene_Name;
						otherwise
							error();
					end
					[SubNet_Score, SN_sid] = sort(SubNet_Score, 'Descend');
					SN_sid(SubNet_Score==0) = [];
					if numel(SN_sid)>MAX_N_SN
						SN_sid = SN_sid(1:MAX_N_SN);
					end
					n_SN = min([numel(res_data.SubNet_List) numel(SN_sid)]);
					Gene_Score = rand(n_gene, 1)*eps; % To make sure stability is not due to selection of all gene all the time
					for gi=1:n_SN
						g_ind = getGeneIndex(res_data.Gene_Name(res_data.SubNet_List{SN_sid(gi)}));
						Gene_Score(g_ind) = Gene_Score(g_ind) + MAX_N_SN - gi + 1;
					end
					[TopGene_Val, TopGene_Ind] = sort(Gene_Score, 'Descend');
					TopGene_Ind(TopGene_Val<1) = 0;
					Part_Marker(step, :) = [mi ni si ri round(res_data.te_auc*100, 3) TopGene_Ind(1:n_TopMarker)'];
					Part_AUC(si, ri) = Part_Marker(step,5);
					step = step + 1;
				end
			end
			save(sav_name, 'Part_Marker', 'Part_AUC', 'Gene_Map');
		end
		Result_Marker = [Result_Marker; Part_Marker];
		Result_AUC(mi,ni,:,:) = Part_AUC;
	end
end
if any(isnan(Result_AUC(:)))
	fprintf('[i] Warning some values are missing.\n');
end
Gene_Name = GE_Data.Gene_Name;

%% Saving
out_name = sprintf('./Collected_NOP_CV%02d_Markers.mat', cv_ind);
save(out_name, 'Result_Marker', 'Result_AUC', 'Gene_Name', 'net_lst', 'method_lst');

% //////////////////// Functions
function g_ind = getGeneIndex(Gene_List)
global Gene_Map
g_ind = zeros(size(Gene_List));
for gi=1:numel(g_ind)
	g_ind(gi) = Gene_Map(Gene_List{gi});
end
end
