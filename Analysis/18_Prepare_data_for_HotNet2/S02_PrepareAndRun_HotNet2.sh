#!/usr/bin/env bash

# Submit: sbatch --job-name=HotNet2 --output=Logs/HN.%J_%a-%N.out --partition=general --qos=long --mem=16GB --time=100:00:00 --ntasks=1 --cpus-per-task=50 S02_PrepareAndRun_HotNet2.sh
# Visualize python2 server.py -i ../../HotNet2_results/<net name>/
set -o xtrace

echo "Begin Date: `date +'%d/%m/%Y - %H:%M:%S'`"

## Initialization
hotnet2=/tudelft.net/staff-bulk/ewi/insy/DBL/aallahyar/tmp/hotnet2_test2/hotnet2
num_cores=30
num_network_permutations=100
num_heat_permutations=1000
#net_name=SyNet
#net_name=STRING
net_name=SyNet-Shuff
mkdir -p ./TMP_Dir/${net_name}
mkdir -p ./HotNet2_results/${net_name}/

# Create network data.
for beta in 0.4 0.5; do
for net_id in ${net_name}_NP003544; do
echo "Creating network data for [$net_id] with [$beta] beta on `date`"
python $hotnet2/makeNetworkFiles.py \
     -e  HotNet2_Input_Files/${net_id}_edgelist.tsv \
     -i  HotNet2_Input_Files/${net_name}_geneindex.tsv \
     -nn ${net_id} \
     -p  ${net_id} \
     -b  $beta \
     -o  TMP_Dir/${net_name} \
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

echo "Running HotNet2 on `date`"
date
python $hotnet2/HotNet2.py \
    -nf  TMP_Dir/${net_name}/${net_name}_NP003544_ppr_0.4.h5 \
         TMP_Dir/${net_name}/${net_name}_NP003544_ppr_0.5.h5 \
    -pnp TMP_Dir/${net_name}/permuted/${net_name}_NP003544_ppr_0.4_##NUM##.h5 \
         TMP_Dir/${net_name}/permuted/${net_name}_NP003544_ppr_0.5_##NUM##.h5 \
    -hf  HotNet2_Input_Files/iCOGS_MD10.0k_heatfile.json \
    -np  $num_network_permutations \
    -hp  $num_heat_permutations \
    -o   HotNet2_results/${net_name}/ \
    -c   $num_cores

echo "End Date: `date +'%d/%m/%Y - %H:%M:%S'`"
echo "Usage statistics for [$SLURM_JOB_ID]:"
sacct -lj $SLURM_JOB_ID
set +o xtrace
