
for policy in RANDOM PUBLICATIONS; do
    for mu in 0.05 0.25 0.5 0.75 0.95; do
        ./scimod-agency mutationParametersData \
            --fprMutationRate=$mu --policy=$policy --nTrials=10
    done
done
