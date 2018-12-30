#!/bin/bash

#SBATCH --job-name=scimod-agency
#SBATCH -t 0-1:00
#SBATCH -n 20
#SBATCH -N 1
#SBATCH --output=logs/scimod-agency_%A_%a.out
#SBATCH --error=logs/scimod-agency_%A_%a.err
#SBATCH -p fast.q

PARAMS=`awk NR==$SLURM_ARRAY_TASK_ID $1`
echo $PARAMS


# ./scimod-agency "/home/scratch/mturner8/Figs3And5" \
#     --fprMutationRate=0.05 \
#     --fprMutationMagnitude=0.01 \
#     --initialFalsePositiveRate=0.05 \
#     --nIter=10000000 \
#     --nTrials=1 \
#     --paramsList=$PARAMS \
#     --syncFPRs=true

# The base parameters for the main paper results.
./scimod-agency "/home/scratch/mturner8/finaldraft-scimod/" \
    --fprMutationRate=0.05 \
    --fprMutationMagnitude=0.01 \
    --initialFalsePositiveRate=0.05 \
    --nIter=10000000 \
    --nTrials=50 \
    --paramsList=$PARAMS

# Same parameters as above except testing baseRate = 0.5.
# ./scimod-agency "/home/scratch/mturner8/scimod-baseRate0.5/" \
#     --fprMutationRate=0.05 \
#     --fprMutationMagnitude=0.01 \
#     --initialFalsePositiveRate=0.05 \
#     --baseRate=0.5 \
#     --nIter=10000000 \
#     --nTrials=50 \
#     --paramsList=$PARAMS


# Base rate = 0.1 with Wright-Fisher selection.
# ./scimod-agency "/home/scratch/mturner8/scimod-wright-fisher/" \
#     --fprMutationRate=0.05 \
#     --fprMutationMagnitude=0.01 \
#     --initialFalsePositiveRate=0.05 \
#     --nIter=10000000 \
#     --selectionMethod=WRIGHT_FISHER \
#     --nTrials=50 \
#     --paramsList=$PARAMS

# Base rate = 0.5 with Wright-Fisher selection.
# ./scimod-agency "/home/scratch/mturner8/scimod-wright-fisher-b=0.5/" \
#     --fprMutationRate=0.05 \
#     --fprMutationMagnitude=0.01 \
#     --initialFalsePositiveRate=0.05 \
#     --nIter=10000000 \
#     --baseRate=0.5 \
#     --selectionMethod=WRIGHT_FISHER \
#     --nTrials=50 \
#     --paramsList=$PARAMS

