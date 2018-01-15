#!/usr/bin/env bash

# Submit: sbatch --job-name=HotNet2 --output=Logs/HN.%J_%a-%N.out --partition=general --qos=long --mem=16GB --time=04:00:00 --ntasks=1 --cpus=36 --cpus-per-task=36 S02_PrepareAndRun_HotNet2.sh SyNet-Shuff
# Visualize python2 ./server.py -i ../../HotNet2_Files/<net name>/Output_DIR/
#set -o xtrace

echo "Begin Date: `date +'%d/%m/%Y - %H:%M:%S'`"

## Initialization
hotnet2=/tudelft.net/staff-bulk/ewi/insy/DBL/aallahyar/tmp/hotnet2_test2/hotnet2
num_cores=36
num_network_permutations=100
num_heat_permutations=1000
src_dir=$1
#SyNet, SyNet-Shuff
rm -rf ./${src_dir}/TMP_DIR/ ${src_dir}/Output_DIR/
mkdir -p ./${src_dir}/TMP_DIR/ ${src_dir}/Output_DIR/

# Create network data.
for beta in 0.3 0.35 0.4 0.45 0.5 0.55 0.6 0.65 0.7; do
echo "Creating network data for [$src_dir] with [$beta] beta on `date`"
python $hotnet2/makeNetworkFiles.py \
     -e  ${src_dir}/Input_DIR/EdgeList.tsv \
     -i  ${src_dir}/Input_DIR/GeneIndex.tsv \
     -nn Network \
     -p  Network \
     -b  $beta \
     -o  ${src_dir}/TMP_DIR/Network-${beta}/ \
     -np $num_network_permutations \
     -c  $num_cores
done

echo "Creating heat from ${src_dir} pvalue file on `date`"
python $hotnet2/makeHeatFile.py \
    scores \
    -hf ${src_dir}/Input_DIR/iCOGS_LogPval_Normalized.tsv \
    -o  ${src_dir}/Input_DIR/iCOGS_LogPval_Normalized.json \
    -n  PVal
    
echo "Creating heat from ${src_dir} #hit/Size file on `date`"
python $hotnet2/makeHeatFile.py \
    scores \
    -hf ${src_dir}/Input_DIR/iCOGS_NHitSize_Normalized.tsv \
    -o  ${src_dir}/Input_DIR/iCOGS_NHitSize_Normalized.json \
    -n  nHitPerSize

echo "Running HotNet2 on `date`"
python $hotnet2/HotNet2.py \
    -nf  ${src_dir}/TMP_DIR/Network-0.3/Network_ppr_0.3.h5 \
         ${src_dir}/TMP_DIR/Network-0.35/Network_ppr_0.35.h5 \
         ${src_dir}/TMP_DIR/Network-0.4/Network_ppr_0.4.h5 \
         ${src_dir}/TMP_DIR/Network-0.45/Network_ppr_0.45.h5 \
         ${src_dir}/TMP_DIR/Network-0.5/Network_ppr_0.5.h5 \
         ${src_dir}/TMP_DIR/Network-0.55/Network_ppr_0.55.h5 \
         ${src_dir}/TMP_DIR/Network-0.6/Network_ppr_0.6.h5 \
         ${src_dir}/TMP_DIR/Network-0.65/Network_ppr_0.65.h5 \
         ${src_dir}/TMP_DIR/Network-0.7/Network_ppr_0.7.h5 \
    -pnp ${src_dir}/TMP_DIR/Network-0.3/permuted/Network_ppr_0.3_##NUM##.h5 \
         ${src_dir}/TMP_DIR/Network-0.35/permuted/Network_ppr_0.35_##NUM##.h5 \
         ${src_dir}/TMP_DIR/Network-0.4/permuted/Network_ppr_0.4_##NUM##.h5 \
         ${src_dir}/TMP_DIR/Network-0.45/permuted/Network_ppr_0.45_##NUM##.h5 \
         ${src_dir}/TMP_DIR/Network-0.5/permuted/Network_ppr_0.5_##NUM##.h5 \
         ${src_dir}/TMP_DIR/Network-0.55/permuted/Network_ppr_0.55_##NUM##.h5 \
         ${src_dir}/TMP_DIR/Network-0.6/permuted/Network_ppr_0.6_##NUM##.h5 \
         ${src_dir}/TMP_DIR/Network-0.65/permuted/Network_ppr_0.65_##NUM##.h5 \
         ${src_dir}/TMP_DIR/Network-0.7/permuted/Network_ppr_0.7_##NUM##.h5 \
    -hf  ${src_dir}/Input_DIR/iCOGS_LogPval_Normalized.json \
         ${src_dir}/Input_DIR/iCOGS_NHitSize_Normalized.json \
    -np  $num_network_permutations \
    -hp  $num_heat_permutations \
    -o   ${src_dir}/Output_DIR/ \
    -c   $num_cores

echo "End Date: `date +'%d/%m/%Y - %H:%M:%S'`"
echo "Usage statistics for [$SLURM_JOB_ID]:"
sacct -lj $SLURM_JOB_ID
# set +o xtrace
