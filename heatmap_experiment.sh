#!/bin/bash

#SBATCH --job-name=scimod-agency
#SBATCH -t 0-0:05
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --output=logs/scimod-agency_%A_%a.out
#SBATCH --error=logs/scimod-agency_%A_%a.err
#SBATCH -p fast.q

PARAMS=`awk NR==$SLURM_ARRAY_TASK_ID heatmap_params_2.txt`

./scimod-agency "/home/scratch/mturner8/scimod/heatmap-highres" \
    --fprMutationRate=0.025 \
    --fprMutationMagnitude=0.01 \
    --initialFalsePositiveRate=0.05 \
    --nTrials=100 \
    --paramsList=$PARAMS

