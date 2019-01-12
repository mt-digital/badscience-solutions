##
# run_funding_experiment.sh
#
# Using this as a runner. Just modify the values in for loops to change
# parameterizations.
#
# Date: 7 Jan 2019
# Author: Matthew A. Turner <maturner01@gmail.com>
#

# 3 Jan 2019: commented out sections are commented so I can 
# make the ten MIXED where policy param X=0.0.


# for policy in MODIFIED_RANDOM MIXED; do
for policy in MIXED; do
    for awardAmount in 10 35 60 85; do
        for npr_fpd_rate in `seq 0.0 0.1 1.0`; do
            # for publishNegativeResultRate in `seq 0.0 0.1 1.0`; do
            # for falsePositiveDetectionRate in `seq 0.0 0.1 1.0`; do
            for policyParam in `seq 0.0 0.1 1.0`; do
                echo "$policy,$awardAmount,$npr_fpd_rate,$npr_fpd_rate,$policyParam"

            done
        done
    done
done
