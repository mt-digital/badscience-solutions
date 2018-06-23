##
# run_funding_experiment.sh
#
# Using this as a runner. Just modify the values in for loops to change
# parameterizations.
#
# Date: 19 Jun 2018
# Author: Matthew A. Turner <maturner01@gmail.com>
#


for policy in RANDOM PUBLICATIONS FPR; do
    for awardAmount in `seq 1 20` 15 `seq 25 5 150`; do # 40 x 3 = 120 total

        echo "Running policy=$policy and awardAmount=$awardAmount"

        qsub -S /bin/bash -q fast.q -cwd -j y -V -l mem_free=96G -pe smp 20 -N funding -o funding.log -e funding.err funding_trials_qsub.sh $policy $awardAmount

    done
done
