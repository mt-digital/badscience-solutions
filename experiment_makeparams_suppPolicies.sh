##
# run_funding_experiment.sh
#
# Using this as a runner. Just modify the values in for loops to change
# parameterizations.
#
# Date: 1 Jan 2019
# Author: Matthew A. Turner <maturner01@gmail.com>
#


for policy in MODIFIED_RANDOM MIXED; do
    for awardAmount in 10 35 60 85; do
        for npr_fpd_rate in `seq 0.0 0.25 1.0`; do
            # for publishNegativeResultRate in `seq 0.0 0.1 1.0`; do
            # for falsePositiveDetectionRate in `seq 0.0 0.1 1.0`; do
            for policyParam in `seq 0.0 0.2 1.0`; do
                echo "$policy,$awardAmount,$npr_fpd_rate,$np_fpd_rate,$policyParam"

            done
        done
    done
done
