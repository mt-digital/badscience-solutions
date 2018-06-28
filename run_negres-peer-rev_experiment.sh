##
# run_funding_experiment.sh
#
# Using this as a runner. Just modify the values in for loops to change
# parameterizations.
#
# Date: 19 Jun 2018
# Author: Matthew A. Turner <maturner01@gmail.com>
#


    # for awardAmount in `seq 1 20` 15 `seq 25 5 150`; do # 40 x 3 = 120 total
for policy in RANDOM PUBLICATIONS; do
    for awardAmount in 1 `seq 5 5 115` ; do
        for publishNegativeResultRate in 0.5 0.9; do
            for falsePositiveDetectionRate in 0.5 0.9; do
                echo \
                    "policy=$policy; awardAmount=$awardAmount;"\
                    "pubNegResRate=$publishNegativeResultRate;" \ 
                    "falsePositiveDetectionRate=$falsePositiveDetectionRate"

                qsub -S /bin/bash -q fast.q -cwd -j y -V -l mem_free=96G -pe \
                    smp 20 -N negres-peer -o negres-peer.log -e negres-peer.err \
                    negres-peer-rev_trials_qsub.sh $policy $awardAmount \
                    $publishNegativeResultRate $falsePositiveDetectionRate

        done
        done
    done
done
