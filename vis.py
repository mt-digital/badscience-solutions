import warnings
from mpl_toolkits.axes_grid1 import AxesGrid
warnings.simplefilter('ignore')

import matplotlib.pyplot as plt
import numpy as np


def plot_means(experiment_data, params=('FPR', '10', '0.10', '0.20')):

    group = experiment_data[params]

    measures = ['nPublications', 'falsePositiveRate', 'falseDiscoveryRate']

    linestyles = ['-', '--', '-.']

    for idx, var in enumerate(measures):
        arr = group[var][:]
        if var == 'falseDiscoveryRate':
            # XXX Hack because it's not being handled correctly in app.d
            # TODO add unit test to catch this nan in app.d
            arr[np.isnan(arr)] = 0.0
        m = arr.mean(axis=0)

        plt.plot(m, label=var, ls=linestyles[idx])

    plt.legend()
