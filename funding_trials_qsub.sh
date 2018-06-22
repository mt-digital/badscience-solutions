#!/bin/bash

policy=$1
awardAmount=$2

./scimod-agency /scratch/mturner8/scimod/fundingExperiment2 \
    --fprMutationRate=0.1 \
    --fprMutationMagnitude=0.01 \
    --policy=$policy \
    --awardAmount=$awardAmount \
    --initialFalsePositiveRate=0.5 \
    --nTrials=100
