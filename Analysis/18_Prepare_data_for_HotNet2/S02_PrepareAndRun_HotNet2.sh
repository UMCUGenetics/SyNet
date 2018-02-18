#!/usr/bin/env bash

# Submit: sbatch --job-name=HotNet2 --output=Logs/HN.%J_%a-%N.out --partition=general --qos=long --mem=16GB --time=04:00:00 --ntasks=1 --cpus=36 --cpus-per-task=36 S02_PrepareAndRun_HotNet2.sh SyNet-Shuff
# Visualize python2 ./server.py -i ../../HotNet2_Files/Output_DIR/<net name>
#set -o xtrace

echo "Begin Date: `date +'%d/%m/%Y - %H:%M:%S'`"

## Initialization
hotnet2=/tudelft.net/staff-bulk/ewi/insy/DBL/aallahyar/tmp/hotnet2_test2/hotnet2
num_cores=36
num_network_permutations=10
num_heat_permutations=1000
src_dir=$1
#SyNet, SyNet-Shuff
rm -rf ./TMP_DIR/${src_dir}/ ./Output_DIR/${src_dir}/
mkdir -p ./TMP_DIR/${src_dir}/ ./Output_DIR/${src_dir}/

# Create network data.
pos_net_lst=""
neg_net_lst=""
for beta in 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7 0.75 0.8; do
echo "Creating network data for [$src_dir] with [$beta] beta on `date`"
python $hotnet2/makeNetworkFiles.py \
     -e  ./Input_DIR/${src_dir}/EdgeList.tsv \
     -i  ./Input_DIR/${src_dir}/GeneIndex.tsv \
     -nn Network \
     -p  Network \
     -b  $beta \
     -o  ./TMP_DIR/${src_dir}/Network-${beta}/ \
     -np $num_network_permutations \
     -c  $num_cores
pos_net_lst=${pos_net_lst}" "./TMP_DIR/${src_dir}/Network-${beta}/Network_ppr_${beta}.h5
neg_net_lst=${neg_net_lst}" "./TMP_DIR/${src_dir}/Network-${beta}/permuted/Network_ppr_${beta}_##NUM##.h5
done
echo -e "\n\n====\nFinal positive network set is:\n"${pos_net_lst}
echo -e "Final negative network set is:\n"${neg_net_lst}"\n\n"

echo "Creating heat from ${src_dir} pvalue file on `date`"
python $hotnet2/makeHeatFile.py \
    scores \
    -hf ./Input_DIR/${src_dir}/iCOGS_LogPval_Normalized.tsv \
    -o  ./Input_DIR/${src_dir}/iCOGS_LogPval_Normalized.json \
    -n  PVal
    
#echo "Creating heat from ${src_dir} #hit/Size file on `date`"
#python $hotnet2/makeHeatFile.py \
#    scores \
#    -hf ./Input_DIR/${src_dir}/iCOGS_NHitSize_Normalized.tsv \
#    -o  ./Input_DIR/${src_dir}/iCOGS_NHitSize_Normalized.json \
#    -n  nHitPerSize

echo "Running HotNet2 on `date`"
python $hotnet2/HotNet2.py \
    -nf  ${pos_net_lst} \
    -pnp ${neg_net_lst} \
    -hf  ./Input_DIR/${src_dir}/iCOGS_LogPval_Normalized.json \
    -np  $num_network_permutations \
    -hp  $num_heat_permutations \
    -o   ./Output_DIR/${src_dir}/ \
    -c   $num_cores
# ./Input_DIR/${src_dir}/iCOGS_NHitSize_Normalized.json \

echo "End Date: `date +'%d/%m/%Y - %H:%M:%S'`"
echo "Usage statistics for [$SLURM_JOB_ID]:"
sacct -lj $SLURM_JOB_ID
# set +o xtrace
