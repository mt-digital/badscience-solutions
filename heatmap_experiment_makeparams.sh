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
    for awardAmount in 10 40 70 100; do
        for publishNegativeResultRate in 0.1 0.25 0.75 0.9; do
            for falsePositiveDetectionRate in 0.1 0.25 0.75 0.9; do
                echo \
                    "$policy,$awardAmount,$publishNegativeResultRate,$falsePositiveDetectionRate" >> heatmap_params.txt
                    
            done
        done
    done
done
