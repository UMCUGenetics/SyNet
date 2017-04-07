function path_str = getPath(data_name)
switch data_name
	case 'ACES'
		path_str = '../../../Dataset/ACES/GE.mat';
	case 'GBM'
		path_str = '../../../Dataset/Cancer_Genomics_Browser/CGB_GBM.mat';
	case 'GBMFl'
		path_str = '../95_Filter_GBM/CGB_GBM_FreqFiltered.mat';
	case 'META'
		path_str = '../94_Pairwise_Evaluation_of_Genes_METABRIC/METABRIC_Filtered.mat';
	case 'HAIBE'
		path_str = '../96_Pairwise_Evaluation_of_Genes_Haibe/Haibe_Filtered.mat';
	case 'CMB'
		path_str = '../108_Evaluation_of_NOPS_On_Networks/CMB_ACES-META-HAIBE.mat';
	case 'SyNet'
		path_str = '../../Gene_Expression_Datasets/SyNet/SyNet_BatchCorrected.mat';
end
end