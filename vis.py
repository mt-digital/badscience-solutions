import warnings
from mpl_toolkits.axes_grid1 import AxesGrid
warnings.simplefilter('ignore')

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def heatmap(experiment_data, policy='FPR', award_amount='10', ax=None,
            measure='falseDiscoveryRate'):

    fpdrs = experiment_data.fpdrs
    pubneg_rates = experiment_data.pubneg_rates

    def order_key(a):
        return float(a)

    fpdrs.sort(key=order_key)
    pubneg_rates.sort(key=order_key)

    by_fpdr_pubneg = experiment_data[policy, award_amount]
    hm_data = np.zeros((len(fpdrs), len(pubneg_rates)))

    for ii, fpdr in enumerate(fpdrs):
        for jj, pubneg_rate in enumerate(pubneg_rates):

            if measure == 'falseDiscoveryRate':
                _data = by_fpdr_pubneg['{}/{}'.format(fpdr, pubneg_rate)
                                       ][measure][:]
                _data[np.isnan(_data)] = 0.0
                _data = _data.mean(axis=0)[-1]

            else:
                _data = by_fpdr_pubneg['{}/{}'.format(fpdr, pubneg_rate)
                                       ][measure][:].mean(axis=0)[-1]

            hm_data[ii, jj] = _data

    xtl = [
        lab if idx in [0, 5, 10] else ''
        for idx, lab in enumerate(fpdrs)
    ]
    ytl = [
        lab if idx in [0, 5, 10] else ''
        for idx, lab in enumerate(pubneg_rates)
    ]
    ax = sns.heatmap(hm_data, vmin=0.0, vmax=1.0, ax=ax,
                     xticklabels=xtl, yticklabels=ytl)
    ax.set_xlabel('false positive detection rate', size=12)
    ax.set_ylabel('negative result pub rate', size=12)
    ax.invert_yaxis()


def plot_means(experiment_data, params=('FPR', '10', '0.10', '0.20')):

    group = experiment_data[params]

    # measures = ['nPublications', 'falsePositiveRate', 'falseDiscoveryRate']
    measures = ['falsePositiveRate', 'falseDiscoveryRate']

    linestyles = ['-', '--', '-.']

    for idx, var in enumerate(measures):
        arr = group[var][:]
        if var == 'falseDiscoveryRate':
            # XXX Hack because it's not being handled correctly in app.d
            # TODO add unit test to catch this nan in app.d
            arr[np.isnan(arr)] = 0.0
        m = arr.mean(axis=0)

        plt.plot(m, label=var, ls=linestyles[idx])

    plt.ylim(0, 1.0)
    plt.legend()
