#!/bin/bash

#SBATCH --job-name=scimod-agency
#SBATCH -t 0-1:00
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --output=logs/scimod-agency_%A_%a.out
#SBATCH --error=logs/scimod-agency_%A_%a.err
#SBATCH -p fast.q

PARAMS=`awk NR==$SLURM_ARRAY_TASK_ID $1`

./scimod-agency "/home/scratch/mturner8/scimod/10M_funds_pubs/" \
    --fprMutationRate=0.025 \
    --fprMutationMagnitude=0.01 \
    --initialFalsePositiveRate=0.05 \
    --nIter=10000000 \
    --nTrials=50 \
    --paramsList=$PARAMS

