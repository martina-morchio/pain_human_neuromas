#!/bin/bash
#SBATCH --partition=gpu
#SBATCH --qos=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=4
#SBATCH --gres=gpu:4

# shell script for submission of Cellbender analysis to HPC cluster
# repeated for each sample individually
# view supplementary materials for parameters used for each sample


cd /mnt/parscratch/users/mdq19mm/nucseq/combined_runs/cellbender
module load Anaconda3/2022.10
source activate /mnt/parscratch/users/mdq19mm/miniconda3/envs/CellBender
module load CUDA/11.7.0

cellbender remove-background \
--input /mnt/parscratch/users/mdq19mm/nucseq/combined_runs/raw_data/N1/count/sample_raw_feature_bc_matrix.h5 \
--output N1_cellbender_feature_bc_matrix.h5 \
--expected-cells 10000 \
--cuda \
--total-droplets-included 50000 \
--fpr 0.01 \
--epochs 150
