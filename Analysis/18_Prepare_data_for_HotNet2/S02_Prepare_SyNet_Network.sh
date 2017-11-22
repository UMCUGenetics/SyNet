#!/usr/bin/env bash

# Submit: sbatch --job-name=HotNet2 --output=Logs/HN.%J_%a-%N.out --partition=general --qos=long --mem=16GB --time=40:00:00 --ntasks=1 --cpus-per-task=32 S02_Prepare_SyNet_Network.sh
# Visualize python2 server.py -i ../../HotNet2_results/
set -o xtrace

echo "Begin Date: `date +'%d/%m/%Y - %H:%M:%S'`"

hotnet2=/tudelft.net/staff-bulk/ewi/insy/DBL/aallahyar/tmp/hotnet2_test2/hotnet2/
num_cores=32
num_network_permutations=5
num_heat_permutations=50

# Create network data.
for beta in 0.1 0.3 0.5 0.7 0.9; do
for net_name in SyNet_NP010000 SyNet_NP020000 SyNet_NP050000 SyNet_NP100000; do
echo "Creating network data for [$net_name] with [$beta] beta on `date`"
echo python $hotnet2/makeNetworkFiles.py \
    -e  HotNet2_Input_Files/${net_name}_edgelist.tsv \
    -i  HotNet2_Input_Files/SyNet_geneindex.tsv \
    -nn ${net_name} \
    -p  ${net_name} \
    -b  $beta \
    -o  Output/SyNet \
    -np $num_network_permutations \
    -c  $num_cores
done
done


echo "Creating heat from iCOGS-10k on `date`"
python $hotnet2/makeHeatFile.py \
    scores \
    -hf HotNet2_Input_Files/iCOGS_MD10.0k_heatfile.tsv \
    -o  HotNet2_Input_Files/iCOGS_MD10.0k_heatfile.json \
    -n  iCOGS_10k
echo "Creating heat from iCOGS-50k on `date`"
python $hotnet2/makeHeatFile.py \
    scores \
    -hf HotNet2_Input_Files/iCOGS_MD50.0k_heatfile.tsv \
    -o  HotNet2_Input_Files/iCOGS_MD50.0k_heatfile.json \
    -n  iCOGS_50k


echo "Running HotNet2 on `date`"
date
python $hotnet2/HotNet2.py \
    -nf  Output/SyNet/SyNet_NP010000_ppr_0.1.h5 \
         Output/SyNet/SyNet_NP010000_ppr_0.3.h5 \
         Output/SyNet/SyNet_NP010000_ppr_0.5.h5 \
         Output/SyNet/SyNet_NP010000_ppr_0.7.h5 \
         Output/SyNet/SyNet_NP010000_ppr_0.9.h5 \
         Output/SyNet/SyNet_NP020000_ppr_0.1.h5 \
         Output/SyNet/SyNet_NP020000_ppr_0.3.h5 \
         Output/SyNet/SyNet_NP020000_ppr_0.5.h5 \
         Output/SyNet/SyNet_NP020000_ppr_0.7.h5 \
         Output/SyNet/SyNet_NP020000_ppr_0.9.h5 \
         Output/SyNet/SyNet_NP050000_ppr_0.1.h5 \
         Output/SyNet/SyNet_NP050000_ppr_0.3.h5 \
         Output/SyNet/SyNet_NP050000_ppr_0.5.h5 \
         Output/SyNet/SyNet_NP050000_ppr_0.7.h5 \
         Output/SyNet/SyNet_NP050000_ppr_0.9.h5 \
         Output/SyNet/SyNet_NP100000_ppr_0.1.h5 \
         Output/SyNet/SyNet_NP100000_ppr_0.3.h5 \
         Output/SyNet/SyNet_NP100000_ppr_0.5.h5 \
         Output/SyNet/SyNet_NP100000_ppr_0.7.h5 \
         Output/SyNet/SyNet_NP100000_ppr_0.9.h5 \
    -pnp Output/SyNet/permuted/SyNet_NP010000_ppr_0.1_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP010000_ppr_0.3_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP010000_ppr_0.5_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP010000_ppr_0.7_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP010000_ppr_0.9_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP020000_ppr_0.1_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP020000_ppr_0.3_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP020000_ppr_0.5_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP020000_ppr_0.7_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP020000_ppr_0.9_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP050000_ppr_0.1_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP050000_ppr_0.3_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP050000_ppr_0.5_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP050000_ppr_0.7_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP050000_ppr_0.9_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP100000_ppr_0.1_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP100000_ppr_0.3_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP100000_ppr_0.5_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP100000_ppr_0.7_##NUM##.h5 \
         Output/SyNet/permuted/SyNet_NP100000_ppr_0.9_##NUM##.h5 \
    -hf  HotNet2_Input_Files/iCOGS_MD10.0k_heatfile.json \
         HotNet2_Input_Files/iCOGS_MD50.0k_heatfile.json \
    -np  $num_network_permutations \
    -hp  $num_heat_permutations \
    -o   HotNet2_results \
    -c   $num_cores

echo "End Date: `date +'%d/%m/%Y - %H:%M:%S'`"
echo "Usage statistics for [$SLURM_JOB_ID]:"
sacct -lj $SLURM_JOB_ID
set +o xtrace
