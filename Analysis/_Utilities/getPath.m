function path_str = getPath(data_name)
switch data_name
	case 'ACES'
		path_str = '../../../../../Dataset/ACES/GE.mat';
	case 'TCGA'
		path_str = '../../Gene_Expression_Datasets/TCGA/TCGA_Combined.mat';
	case 'UTCGA'
		path_str = '../../Gene_Expression_Datasets/TCGA/TCGA_Unnormalized.mat';
% 	case 'GBM'
% 		path_str = '../../../Dataset/Cancer_Genomics_Browser/CGB_GBM.mat';
% 	case 'GBMFl'
% 		path_str = '../95_Filter_GBM/CGB_GBM_FreqFiltered.mat';
% 	case 'META'
% 		path_str = '../94_Pairwise_Evaluation_of_Genes_METABRIC/METABRIC_Filtered.mat';
% 	case 'HAIBE'
% 		path_str = '../96_Pairwise_Evaluation_of_Genes_Haibe/Haibe_Filtered.mat';
% 	case 'CMB'
% 		path_str = '../108_Evaluation_of_NOPS_On_Networks/CMB_ACES-META-HAIBE.mat';
	case 'SyNet'
		path_str = '../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat';

	%% Biological networks
	case 'STRING'
		path_str = '../../Networks/STRING/STRING_combined_score_GN_WithWeights.txt';
	case 'STR.CoExpr'
		path_str = '../../Networks/STRING/STRING_coexpression_GN_WithWeights.txt';
	case 'STR.CoOccr'
		path_str = '../../Networks/STRING/STRING_cooccurence_GN_WithWeights.txt';
	case 'STR.DB'
		path_str = '../../Networks/STRING/STRING_database_GN_WithWeights.txt';
	case 'STR.Exprm'
		path_str = '../../Networks/STRING/STRING_experimental_GN_WithWeights.txt';
	case 'STR.Fus'
		path_str = '../../Networks/STRING/STRING_fusion_GN_WithWeights.txt';
	case 'STR.Neigh'
		path_str = '../../Networks/STRING/STRING_neighborhood_GN_WithWeights.txt';
	case 'STR.TxtMn'
		path_str = '../../Networks/STRING/STRING_textmining_GN_WithWeights.txt';
	case 'KEGG'
		path_str = '../../Networks/KEGG/KEGG.txt';
	case 'MSigDB'
		path_str = '../../Networks/MSigDB/MSigDB.tsv';
	case 'HPRD'
		path_str = '../../Networks/HPRD/HPRD.tsv';
	case 'I2D'
		path_str = '../../Networks/I2D/I2D.tsv';
    case 'HBEpith'
		path_str = '../../Networks/HumanBase/HumanBase_mammary_epithelium_Top1Mil.tsv';
    case 'HBGland'
		path_str = '../../Networks/HumanBase/HumanBase_mammary_gland_Top1Mil.tsv';
	otherwise
		error();
end
end