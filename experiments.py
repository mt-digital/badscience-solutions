from subprocess import call


def mutation_parameters():

    fpr_mutation_rates = [0.5]  # [0.1, 0.25, 0.5, 0.75, 0.9]
    fpr_mutation_magnitudes = [0.01]  # , 0.05, 0.1, 0.25, 0.5]

    for rate in fpr_mutation_rates:
        for mag in fpr_mutation_magnitudes:
            call(
                "./scimod-agency mutationParametersData "
                "--fprMutationRate={} --fprMutationMagnitude={}".format(
                    rate, mag
                 ),
                shell=True
             )
