#!/bin/bash

#SBATCH --job-name=scimod-agency
#SBATCH -t 0-2:30
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --output=logs/json_to_hdf.out
#SBATCH --error=logs/json_to_hdf.err
#SBATCH -p fast.q

python json_to_hdf.py $1 $2

