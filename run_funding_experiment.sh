##
# run_funding_experiment.sh
#
# 
#
# Date: 19 Jun 2018
# Author: Matthew A. Turner <maturner01@gmail.com>
#
# Example:
#
# qsub -S /bin/bash -q fast.q -cwd -j y -V -l mem_free=96G \
#     -pe smp 20 -N mutation_params -o params.log -e params.err \
#     ./run_funding_experiment.sh testData


for policy in RANDOM PUBLICATIONS FPR; do
    for awardAmount in 10 50 100; do

        echo "Running policy=$policy and awardAmount=$awardAmount"

        ./scimod-agency $1 \
            --fprMutationRate=0.25 \
            --fprMutationMagnitude=0.01 \
            --policy=$policy \
            --awardAmount=$awardAmount \
            --nTrials=3
    done
done
