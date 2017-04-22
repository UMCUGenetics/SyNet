clc;
clear;

%% Initialization
result_path = './Results_Files/';
net_lst = {'Corr', 'AbsCorr', 'Random', 'String', 'Multinet', 'KEGG', 'MSigDB'
	'DSN-SyNet-T00500'
	'DSN-SyNet-T01000'
	'DSN-SyNet-P99.999'
	'DSN-SyNet-T05000'
	'DSN-SyNet-T10000'
	'DSN-SyNet-P99.99'
	'DSN-SyNet-T20000'
	'DSN-SyNet-P99.90'
	};
n_net = numel(net_lst);
method_lst = {'iPark', 'iChuang', 'iTaylor', 'Feral', 'Single'};
n_met = numel(method_lst);
cv_id = '170412';
n_study = 14;
n_rep = 20;
MAX_N_SUBNET = 500;
n_Top = 100;

%% Generate gene map
load('../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat', 'Gene_Name');
n_gene = numel(Gene_Name);
global Gene_Map
Gene_Map = containers.Map(Gene_Name, 1:n_gene);
n_row = prod([n_met n_net n_study n_rep]);
Result_Info = zeros(n_row, 5 + n_Top); % Method_Index, Network_Index, Study_Index, Repeat_Index, Test_AUC
Result_AUC = nan(n_met, n_net, n_study, n_rep);
step = 1;
for si=1:n_study
	for ri=1:n_rep
		for ni=1:n_net
			for mi=1:n_met
				net_name = strrep(net_lst{ni}, 'SyNet', sprintf('SyNetS%02d',si));
				res_name = [result_path 'CID-' cv_id 'Si' num2str(si,'%02d') 'Ri' num2str(ri,'%02d') '*_NET-' net_name '_*MSN-' num2str(MAX_N_SUBNET) '_MTN-' method_lst{mi} '.mat'];
				file_info = dir(res_name);
				if numel(file_info)~=1
					fprintf('*** Ignoring [%s]\n', res_name);
					continue;
				end
				fprintf('Reading from: [%s]\n', file_info.name);
				res_data = load([result_path file_info.name]);
				Gene_Score = rand(n_gene, 1)*eps; % To make sure stability is not due to selection of all gene all the time
				feat_scr = abs(res_data.B(:,res_data.fit.IndexMinMSE));
				[feat_val, feat_sid] = sort(feat_scr, 'Descend');
				feat_sid(feat_val==0) = [];
				if numel(feat_sid)>MAX_N_SUBNET
					feat_sid = feat_sid(1:MAX_N_SUBNET);
				end
				n_feat = numel(feat_sid);
				switch method_lst{mi}
					case {'iPark' 'iChuang' 'iTaylor' 'Feral'}
						for fi=1:n_feat
							g_ind = getGeneIndex(res_data.Gene_Name(res_data.SubNet_List{feat_sid(fi)}));
							Gene_Score(g_ind) = Gene_Score(g_ind) + MAX_N_SUBNET - fi + 1;
						end
					case 'Single'
						g_ind = getGeneIndex(res_data.Gene_Name(feat_sid));
						Gene_Score(g_ind) = MAX_N_SUBNET:-1:MAX_N_SUBNET-n_feat+1;
					otherwise
						error();
				end
				[TopGene_Val, TopGene_Ind] = sort(Gene_Score, 'Descend');
				TopGene_Ind(TopGene_Val<1) = 0;
				Result_Info(step, :) = [mi ni si ri round(res_data.te_auc*100, 3) TopGene_Ind(1:n_Top)'];
				Result_AUC(mi,ni,si,ri) = Result_Info(step,5);
				step = step + 1;
			end
		end
	end
end
if any(isnan(Result_AUC(:)))
	fprintf('[i] Warning some values are missing.\n');
end
save('./Collected_NOP_results.mat', 'Result_Info', 'Result_AUC', 'Gene_Name', 'net_lst', 'method_lst');

% //////////////////// Functions 
function g_ind = getGeneIndex(Gene_List)
global Gene_Map

g_ind = zeros(size(Gene_List));
for gi=1:numel(g_ind)
	g_ind(gi) = Gene_Map(Gene_List{gi});
end
end
