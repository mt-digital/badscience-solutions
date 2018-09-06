#!/bin/bash

#SBATCH --job-name=scimod-agency
#SBATCH -t 0-2:30
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --output=logs/scimod-agency_%A_%a.out
#SBATCH --error=logs/scimod-agency_%A_%a.err
#SBATCH -p fast.q

python json_to_hdf.py $1 $2

