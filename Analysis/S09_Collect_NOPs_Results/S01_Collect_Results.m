clc;
clear;

%% Initialization
result_path = '../S08_Perform_NOPs_Restricted_Network/Results_Files/';
net_lst = {
	'SyNet-T10000', 'MinSyn-T10000', 'AvgSyn-T10000', 'CrSyn-T10000', 'CrMinSyn-T10000' ...
	'Corr-T10000','AbsCorr-T10000','Random-T10000', ...
	};
n_net = numel(net_lst);
method_lst = {'TNMC' 'TLEx' 'Lasso' 'iPark' 'RI-iPark' 'Feral'  'RI-Feral'}; % 'Regress' 'RegAG' 'LExAG'
n_met = numel(method_lst);
n_study = 12;
n_rep = 10;
MAX_N_SN = 500;
n_TopMarker = 500;

%% Generate gene map
GE_Data = load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name', 'Study_Name');
n_gene = numel(GE_Data.Gene_Name);
global Gene_Map
Gene_Map = containers.Map(GE_Data.Gene_Name, 1:n_gene);
n_row = prod([n_met n_net n_study n_rep]);
Result_Info = zeros(n_row, 5 + n_TopMarker);
% 1:Method_Index, 2:Network_Index, 3:Study_Index, 4:Repeat_Index, 5:Test_AUC
Result_AUC = nan(n_met, n_net, n_study, n_rep);
step = 1;
for ni=1:n_net
	for mi=1:n_met
		for si=1:n_study
			for ri=1:n_rep
				net_name = net_lst{ni};
				%net_name = strrep(net_lst{ni}, 'SyNet', sprintf('SyNetS%02d',si));
				res_name = [result_path 'DID_Si' num2str(si,'%02d') '-Ri' num2str(ri,'%03d') '_' net_name '_*MSN-' num2str(MAX_N_SN,'%03d') '_MTN-' method_lst{mi} '.mat'];
				file_info = dir(res_name);
				if numel(file_info)~=1
					fprintf('*** Ignoring [%s]\n', res_name);
					continue;
				end
				fprintf('Reading from: [%s]\n', file_info.name);
				res_data = load([result_path file_info.name]);
				switch method_lst{mi}
					case 'TNMC'
						SubNet_Score = n_gene:-1:1;
					case 'Lasso'
						SubNet_Score = abs(res_data.B(:,res_data.fit.IndexMinMSE));
						res_data.SubNet_List = num2cell(1:numel(res_data.Gene_Name))';
					case 'LExAG'
						SubNet_Score = abs(res_data.B(:,res_data.fit.IndexMinMSE));
						res_data.SubNet_List = num2cell(1:n_gene)';
						res_data.Gene_Name = GE_Data.Gene_Name;
					otherwise
						SubNet_Score = res_data.SubNet_Score;
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
				Result_Info(step, :) = [mi ni si ri round(res_data.te_auc*100, 3) TopGene_Ind(1:n_TopMarker)'];
				Result_AUC(mi,ni,si,ri) = Result_Info(step,5);
				step = step + 1;
			end
		end
	end
end
if any(isnan(Result_AUC(:)))
	fprintf('[i] Warning some values are missing.\n');
end
Gene_Name = GE_Data.Gene_Name;
save('./Collected_NOP_results.mat', 'Result_Info', 'Result_AUC', 'Gene_Name', 'net_lst', 'method_lst');

% //////////////////// Functions
function g_ind = getGeneIndex(Gene_List)
global Gene_Map

g_ind = zeros(size(Gene_List));
for gi=1:numel(g_ind)
	g_ind(gi) = Gene_Map(Gene_List{gi});
end
end
