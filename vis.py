import warnings
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.backends.backend_pdf import PdfPages
warnings.simplefilter('ignore')

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


current_palette = sns.color_palette('viridis')
sns.set_palette(current_palette)


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

            try:
                if measure == 'falseDiscoveryRate':
                        _data = by_fpdr_pubneg['{}/{}'.format(fpdr, pubneg_rate)
                                               ][measure][:]
                        _data[np.isnan(_data)] = 0.0
                        _data = _data.mean(axis=0)[-1]

                else:
                    _data = by_fpdr_pubneg['{}/{}'.format(fpdr, pubneg_rate)
                                           ][measure][:].mean(axis=0)[-1]
            except:
                print('fpdr: {}, pubneg_rate: {}'.format(fpdr, pubneg_rate))
                _data = np.nan

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
                     xticklabels=xtl, yticklabels=ytl, square=True,
                     cmap='viridis')
    ax.set_xlabel('false positive detection rate', size=12)
    ax.set_ylabel('negative result pub rate', size=12)
    ax.invert_yaxis()

    ax.set_title('Policy: {}, $G={}$'.format(policy, award_amount))


def plot_means(experiment_data, params=('FPR', '10', '0.10', '0.20'),
               ax=None, equal=True):

    group = experiment_data[params]

    # measures = ['nPublications', 'falsePositiveRate', 'falseDiscoveryRate']
    measures = ['falsePositiveRate', 'falseDiscoveryRate']

    linestyles = ['-', '--', '-.']

    if ax is None:
        fig, ax = plt.subplots()

    lines = [None] * 3
    for idx, var in enumerate(measures):
        arr = group[var][:]
        if var == 'falseDiscoveryRate':
            # XXX Hack because it's not being handled correctly in app.d
            # TODO add unit test to catch this nan in app.d
            arr[np.isnan(arr)] = 0.0
        m = arr.mean(axis=0)

        lines[idx] = ax.plot(m, label=var, ls=linestyles[idx])

    ax.set_ylim(0, 1.0)
    ax.set_ylabel('false positive/discovery rate')

    npub_arr = group['nPublications'][:].mean(axis=0)
    ax2 = ax.twinx()
    lines[2] = ax2.plot(
        npub_arr, label='nPublications', ls=linestyles[-1], color='g'
    )
    ax2.set_ylabel('nPublications', color='g')
    ax2.tick_params('y', colors='g')

    # labels = [l[0].get_label() for l in lines]
    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    labels = ['FPR', 'FDR', 'nPubs']
    ax2.legend(lines + lines2, labels, loc=(0.5, 0.65), fontsize=6)

    ax.set_title(
        'Policy={}\nG={}\nfpdr={}\npubneg_rate={}'.format(*params), fontsize=8
    )


def heatmaps_and_convergence_check(experiment_data):

    award_amounts = [int(amt) for amt in experiment_data.award_amounts]
    award_amounts.sort()
    with PdfPages('heatmaps_and_conv.pdf') as pdf:
        # for award_amount in [award_amounts[0]]:
        for award_amount in award_amounts:
            fig, axes = plt.subplots(3, sharex=True, figsize=(8.5, 11))
            for idx, policy in enumerate(experiment_data.policies):
                heatmap(experiment_data, policy, award_amount, axes[idx])
            pdf.savefig(fig)

        fpdrs_negrates = ['0.20', '0.80', '0.90']
        for award_amount in award_amounts:
            for idx, policy in enumerate(experiment_data.policies):
                fig, axes = plt.subplots(3, 3, figsize=(8.5, 11))
                for fpdr_idx, fpdr in enumerate(fpdrs_negrates):
                    for pubneg_idx, pubneg_rate in enumerate(fpdrs_negrates):
                        try:
                            plot_means(
                                experiment_data,
                                params=(policy, award_amount, fpdr, pubneg_rate),
                                ax=axes[fpdr_idx, pubneg_idx]
                            )
                        # In case a parameter setting got dropped, don't plot.
                        except KeyError:
                            pass
                pdf.savefig(fig)
