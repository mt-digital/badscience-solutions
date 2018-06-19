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
#     -pe smp 20 -N mutation_funding -o funding.log -e funding.err \
#     ./run_funding_experiment.sh testData


for policy in RANDOM PUBLICATIONS FPR; do
    for awardAmount in  `seq 10 5 205`; do # 40 x 3 = 120 total

        echo "Running policy=$policy and awardAmount=$awardAmount"

        qsub -S /bin/bash -q fast.q -cwd -j y -V -l mem_free=96G -pe smp 20 -N funding -o funding.log -e funding.err funding_trials_qsub.sh $policy $awardAmount

    done
done
