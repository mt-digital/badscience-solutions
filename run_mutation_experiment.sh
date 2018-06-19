##
# run_mutation_experiment.sh
#
# Run parameterization experiments of mutation 
#
# Date: 18 Jun 2018
# Author: Matthew A. Turner <maturner01@gmail.com>
#
# Example:
#
# qsub -S /bin/bash -q fast.q -cwd -j y -V -l mem_free=96G \
#     -pe smp 20 -N mutation_params -o params.log -e params.err \
#     ./run_mutation_experiment.sh testData


for policy in RANDOM PUBLICATIONS; do
    for mu_mag in 0.01 0.05 0.1 0.15 0.2 0.25; do
        for mu in 0.01 0.05 0.25 0.5 0.75 0.95; do

            echo "Running policy=$policy; mu_alpha=$mu; mu_mag=$mu_mag"

            ./scimod-agency $1 \
                --fprMutationRate=$mu \
                --policy=$policy \
                --fprMutationMagnitude=$mu_mag \
                --nTrials=20
        done
    done
done
