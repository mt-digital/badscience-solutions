#!/bin/bash
#
# Wrapper around scimod-agency binary that Sun Grid Engine will not run bare.
# 
# Date: 19 Jun 2018
# Author: Matthew A. Turner <maturner01@gmail.com>
#

policy=$1
awardAmount=$2

./scimod-agency /scratch/mturner8/scimod/fundingExperiment \
    --fprMutationRate=0.025 \
    --fprMutationMagnitude=0.01 \
    --policy=$policy \
    --awardAmount=$awardAmount \
    --initialFalsePositiveRate=0.05 \
    --nTrials=100
