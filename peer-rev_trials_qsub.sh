#!/bin/bash
#
# Wrapper around scimod-agency binary that Sun Grid Engine will not run bare.
# 
# Date: 19 Jun 2018
# Author: Matthew A. Turner <maturner01@gmail.com>
#

policy=$1
awardAmount=$2
falsePositiveDetectionRate=$3

./scimod-agency /scratch/mturner8/scimod/peer-review-gte0p8/ \
    --fprMutationRate=0.025 \
    --fprMutationMagnitude=0.01 \
    --policy=$policy \
    --awardAmount=$awardAmount \
    --initialFalsePositiveRate=0.05 \
    --publishNegativeResultRate=0.0 \
    --falsePositiveDetectionRate=$falsePositiveDetectionRate \
    --nTrials=100
