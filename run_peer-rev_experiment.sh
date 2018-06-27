##
# run_peer-rev_experiment.sh
#
# Using this as a runner. Just modify the values in for loops to change
# parameterizations.
#
# Date: 27 Jun 2018
# Author: Matthew A. Turner <maturner01@gmail.com>
#


    # for awardAmount in `seq 1 20` 15 `seq 25 5 150`; do # 40 x 3 = 120 total
for policy in RANDOM PUBLICATIONS; do
    for awardAmount in 1 `seq 5 5 115` ; do # 40 x 3 = 120 total
        for falsePositiveDetectionRate in 0.25 0.5 0.75 1.0; do
            echo "policy=$policy; awardAmount=$awardAmount; falsePositiveDetectionRate=$falsePositiveDetectionRate"

            qsub -S /bin/bash -q fast.q -cwd -j y -V -l mem_free=96G -pe \
                smp 20 -N peer -o peer.log -e peer.err \
                peer-rev_trials_qsub.sh $policy $awardAmount $falsePositiveDetectionRate

        done
    done
done
