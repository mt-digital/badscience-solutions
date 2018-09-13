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
for policy in FPR RANDOM PUBLICATIONS; do
    for awardAmount in 10 35 60 85; do
        for publishNegativeResultRate in `seq 0.0 0.1 1.0`; do
            for falsePositiveDetectionRate in `seq 0.0 0.1 1.0`; do
                echo \
                    "$policy,$awardAmount,$publishNegativeResultRate,$falsePositiveDetectionRate"
                    
            done
        done
    done
done
