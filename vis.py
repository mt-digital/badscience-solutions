import warnings

from glob import glob
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.backends.backend_pdf import PdfPages
from itertools import repeat
warnings.simplefilter('ignore')

import json
import matplotlib.pyplot as plt
plt.style.use('seaborn-deep')
import matplotlib.lines as mlines

import matplotlib.pyplot as plt
import numpy as np
import os
import seaborn as sns


# current_palette = sns.color_palette('viridis')
# sns.set_palette(current_palette)


def heatmap(experiment_data, policy='FPR', award_amount='10', ax=None,
            measure='falseDiscoveryRate', ylabel=True, colorbar=False, **kwargs):

    fpdrs = experiment_data.fpdrs
    pubneg_rates = experiment_data.pubneg_rates

    def order_key(a):
        return float(a)

    fpdrs.sort(key=order_key)
    pubneg_rates.sort(key=order_key)

    by_fpdr_pubneg = experiment_data[policy, award_amount]
    hm_data = np.zeros((len(fpdrs), len(pubneg_rates)))

    for ii, pubneg_rate in enumerate(pubneg_rates):
        for jj, fpdr in enumerate(fpdrs):

            try:
                if measure == 'falseDiscoveryRate':
                        _data = by_fpdr_pubneg['{}/{}'.format(pubneg_rate, fpdr)
                                               ][measure][:]
                        _data[np.isnan(_data)] = 0.0
                        _data = _data.mean(axis=0)[-1]

                else:
                    _data = by_fpdr_pubneg['{}/{}'.format(pubneg_rate, fpdr)
                                           ][measure][:].mean(axis=0)[-1]
            except:
                print('fpdr: {}, pubneg_rate: {}'.format(pubneg_rate, fpdr))
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
                     cmap='viridis', cbar=colorbar, **kwargs)

    ax.set_xlabel(r'efficacy of peer review, $r$', size=14)
    if ylabel:
        ax.set_ylabel(r'pub. rate of neg. results, $p$', size=14)

    ax.invert_yaxis()

    ax.set_title('Policy: {}, $G={}$'.format(
        PAPER_POLICY_LOOKUP[policy], award_amount)
    )


def heatmaps(experiment_data, award_amount,
             measure='falseDiscoveryRate', save_path=None, figsize=(10, 6)):

    fig, axes = plt.subplots(1, 3, sharey=True, figsize=figsize)
    cbar_ax = fig.add_axes([.99, .275, .03, .45])

    colorbar = None
    for idx, policy in enumerate(POLICIES):
        ax = axes[idx]
        if idx == 0:
            label = r'$F$' if measure == 'falseDiscoveryRate' else r'$\overline{\alpha}$'
            heatmap(experiment_data, policy, ax=ax, award_amount=award_amount,
                    measure=measure, colorbar=True,
                    cbar_ax=cbar_ax)
        else:
            heatmap(experiment_data, policy, ax=ax, award_amount=award_amount,
                    measure=measure, ylabel=False)
        # if idx == 2:
            # heatmap(experiment_data, policy, ax=ax,
            #         measure=measure, colorbar=True)
    ax.figure.axes[-1].set_ylabel(label, size=18)
    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path, bbox_inches='tight')


def plot_means(experiment_data,
               params=('FPR', '10', '0.10', '0.20'),
               ax=None, equal=True, publication_measure='meanPublications'):

    group = experiment_data[params]

    # measures = ['nPublications', 'falsePositiveRate', 'falseDiscoveryRate']
    measures = ['falsePositiveRate', 'falseDiscoveryRate']

    # linestyles = ['-', '--', '-.']
    linestyles = ['-', '-', '-']

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

    ax.set_ylim(0, 1.05)
    ax.set_ylabel('false positive/discovery rate')

    npub_arr = group[publication_measure][:].mean(axis=0)
    ax2 = ax.twinx()
    lines[2] = ax2.plot(
        npub_arr, label='pubs', ls=linestyles[-1], color='g'
    )
    ax2.set_ylabel('pubs', color='g')
    ax2.tick_params('y', colors='g')

    lines, labels = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    labels = [r'$\overline{\alpha}$', 'F', 'tot. pubs']
    ax2.legend(lines + lines2, labels, loc=(0.1, 0.55), ncol=3, handlelength=1,
               fontsize=9)

    ax.set_xticks([0, 500, 1000])
    ax.set_xticklabels(['0', '5e6', '1e7'])
    ax.set_xlabel('Iteration')

    params = [ PAPER_POLICY_LOOKUP[params[0]] ] + list(params[1:])
    ax.set_title(
        # 'Policy={}\nG={}\n$n$={}\n$d$={}'.format(*params),
        'Policy={}, G={}, $n$={}, $d$={}'.format(*params),
        fontsize=10
    )


def heatmaps_and_convergence_check(experiment_data,
                                   heatmaps=True,
                                   publication_measure='meanPublications',
                                   save_path='heatmaps_and_conv.pdf',
                                   save_dir=None,
                                   figsize=(3, 3)):
    '''
    Create heatmaps and demonstrations of convergence for a selection of
    parameter settings. If save_dir is not provided, a PDF of plots of behavior
    across all three policies, all four grant funding amounts, and three
    values each for fpdr and negative publishing rate. If save_dir is provided,
    a separate PDF of each timeseries is created and saved to save_dir with
    filename convention "policy_FPR-G_10-fpdr_0.2-npr_0.1.pdf".
    '''

    award_amounts = [int(amt) for amt in experiment_data.award_amounts]
    award_amounts.sort()

    fpdrs_negrates = ['0.20', '0.80', '0.90']
    # Without save directory, create a mega-PDF of all timeseries.
    if save_dir == None:
        with PdfPages(save_path) as pdf:
            # for award_amount in [award_amounts[0]]:
            if heatmaps:
                for award_amount in award_amounts:
                    fig, axes = plt.subplots(3, sharex=True, figsize=(8.5, 11))
                    for idx, policy in enumerate(experiment_data.policies):
                        heatmap(experiment_data, policy, award_amount, axes[idx])
                    pdf.savefig(fig)

            for award_amount in award_amounts:
                for idx, policy in enumerate(experiment_data.policies):
                    fig, axes = plt.subplots(3, 3, figsize=(8.5, 11))
                    for pubneg_idx, pubneg_rate in enumerate(fpdrs_negrates):
                        for fpdr_idx, fpdr in enumerate(fpdrs_negrates):
                            try:
                                plot_means(
                                    experiment_data,
                                    params=(policy, award_amount, pubneg_rate, fpdr),
                                    ax=axes[pubneg_idx, fpdr_idx],
                                    publication_measure=publication_measure
                                )
                            # In case a parameter setting got dropped, don't plot.
                            except KeyError:
                                pass
                    pdf.savefig(fig)
    else:
        # No heatmaps, but do create PDF for each policy with name as above.
        import matplotlib as mpl
        mpl.style.use('default')
        for award_amount in award_amounts:
            for idx, policy in enumerate(experiment_data.policies):
                for pubneg_idx, pubneg_rate in enumerate(fpdrs_negrates):
                    for fpdr_idx, fpdr in enumerate(fpdrs_negrates):
                        if pubneg_idx < 1:
                            fig = plt.figure(figsize=figsize)
                            plot_means(
                                experiment_data,
                                params=(policy,
                                        award_amount, pubneg_rate, fpdr),
                                publication_measure=publication_measure,
                                ax=plt.gca()
                            )
                            plt.savefig(
                                os.path.join(
                                    save_dir,
                                    'policy_{}-G_{}-npr_{}-fpdr_{}.pdf'.format(
                                    policy, award_amount, pubneg_rate, fpdr)
                                )
                            )


def policy_diff_heatmap(experiment_data, award_amount,
                        measure='falseDiscoveryRate'):
    '''
    Plot a heatmap of the difference between the publications and random
    policies to see how their outcomes differ.
    '''

    fpdrs = experiment_data.fpdrs
    pubneg_rates = experiment_data.pubneg_rates

    def order_key(a):
        return float(a)

    fpdrs.sort(key=order_key)
    pubneg_rates.sort(key=order_key)

    publications = experiment_data['PUBLICATIONS', award_amount]
    random = experiment_data['RANDOM', award_amount]

    hm_data = np.zeros((len(fpdrs), len(pubneg_rates)))
    for ii, pubneg_rate in enumerate(pubneg_rates):
        for jj, fpdr in enumerate(fpdrs):

            try:
                random_vec = random[
                    '{}/{}'.format(pubneg_rate, fpdr)][measure][:]
                publication_vec = publications[
                    '{}/{}'.format(pubneg_rate, fpdr)
                        ][measure][:]

                # XXX Still have to have hack to get rid of nan, fix this!
                random_vec[np.isnan(_data)] = 0.0
                publication_vec[np.isnan(_data)] = 0.0

                random_mean = random_vec.mean(axis=0)[-1]
                publication_mean = publication_vec.mean(axis=0)[-1]
                _data = random_mean - publication_mean
            except:
                print('pubneg_rate: {}, fpdr: {}'.format(pubneg_rate, fpdr))
                _data = 0.0

            hm_data[ii, jj] = _data

    xtl = [
        lab if idx in [0, 5, 10] else ''
        for idx, lab in enumerate(pubneg_rates)
    ]
    ytl = [
        lab if idx in [0, 5, 10] else ''
        for idx, lab in enumerate(fpdrs)
    ]
    ax = sns.heatmap(hm_data,
                     xticklabels=xtl, yticklabels=ytl, square=True,
                     cmap='seismic', center=0.0, robust=True)
    ax.set_ylabel('Negative result pub rate (NPR)', size=12)
    ax.set_xlabel('False positive detection rate (FPDR)', size=12)
    ax.invert_yaxis()
    measure_acr = 'FDR' if measure == 'falseDiscoveryRate' else 'FPR'
    ax.figure.axes[-1].set_ylabel(
        '${0}_{{RANDOM}} - {0}_{{PUBLICATIONS}}$'.format(measure_acr))

    # ax.set_title('FDR difference between publication and random policies')

    return hm_data

def measure_vs_pubparams(experiment_data,
                         award_amounts=['10', '35', '60', '85'],
                         fpdr_negres_vals=['0.00', '0.10', '0.20', '0.30',
                                           '0.40', '0.50', '0.60', '0.70',
                                           '0.80', '0.90', '1.00'],
                         other_val='0.00',
                         param='NPR',
                         policy='PUBLICATIONS',
                         measure='falseDiscoveryRate',
                         ax=None,
                         xlabel=False,
                         ylabel=False,
                         legend=False,
                         title=False):
    '''
    Plot either FPR or FDR for different award amounts over all publishing
    parameter values; parameter is either false positive discovery rate ('FPDR')
    or negative (result) publishing rate ('NPR').
    '''
    import matplotlib.pyplot as plt
    plt.style.use('seaborn-deep')

    if ax is None:
        fig, ax = plt.subplots()

    for idx, award_amount in enumerate(award_amounts):
        if param == 'NPR':
            data = np.array(
                [experiment_data[policy, award_amount, npr, other_val][measure][:]
                 for npr in fpdr_negres_vals]
            )
        elif param == 'FPDR':
            data = np.array(
                [experiment_data[policy, award_amount, other_val, fpdr][measure][:]
                 for fpdr in fpdr_negres_vals]
            )
        else:
            raise ValueError('param must be NPR or FPDR')

        data[np.isnan(data)] = 0.0

        dm = data.mean(axis=1)[:, -1]
        ax.plot(dm, 'o-', label='$G={}$'.format(award_amount))

    ax.set_xticks([0, 5, 10])
    ax.set_xticklabels(['0.0', '0.5', '1.0'], fontsize=14)

    if xlabel:
        xlab = (
            r'efficacy of peer review, $r$'
            if param == 'FPDR'
            else r'pub. rate of neg. results, $p$'
        )
        ax.set_xlabel(xlab, size=14)

    if 'Publications' not in measure:
        ax.set_ylim(-0.05, 1.05)
        ax.set_yticks(np.arange(0, 1.01, 0.25))
    ax.yaxis.set_tick_params(labelsize=14)
    ax.grid(True, axis='both')

    if ylabel:
        ax.set_ylabel(ylabel, size=18)
    if title:
        policy_in_pub = PAPER_POLICY_LOOKUP[policy]
        ax.set_title('{} award policy'.format(policy_in_pub),
                     size=16)

    if legend:
        ax.legend(fontsize=12)


def many_measure_vs_subparams(experiment_data, param='FPDR',
                              figsize=(10, 5.5), save_path=None):

    fig, axes = plt.subplots(2, 3, sharey=True, sharex=True, figsize=figsize)
    ax_idx = 0
    measures = ['falseDiscoveryRate', 'falsePositiveRate']

    for m_idx, measure in enumerate(measures):
        ylabel = r'$F$' if m_idx == 0 else r'$\overline{\alpha}$'
        xlabel = False
        title = True
        legend = False
        for p_idx, policy in enumerate(POLICIES):
            if ax_idx == 3:
                legend = True
            if p_idx == 1:
                ylabel = False
            if m_idx == 1:
                xlabel = True
                title = False

            measure_vs_pubparams(experiment_data, measure=measure,
                                 param=param, policy=policy,
                                 ax=axes.flatten()[ax_idx],
                                 legend=legend, ylabel=ylabel, xlabel=xlabel,
                                 title=title)
            legend = False
            ax_idx += 1

    plt.tight_layout()
    if save_path is not None:
        plt.savefig(save_path)


POLICIES = ['PUBLICATIONS', 'RANDOM', 'FPR']
PAPER_POLICY_LOOKUP = {
    'PUBLICATIONS': 'PH',
    'RANDOM': 'RA',
    'FPR': 'MI',
    'MIXED': 'MIXED'
}


def policies_timeseries(experiment_data,
                        award_amounts=['10', '35', '60', '85'],
                        lcs=['blue', 'red', 'black'],
                        figsize=(7, 5),
                        save_path=None):

    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True, figsize=figsize)

    for g_idx, award_amount in enumerate(award_amounts):
        ax = axes.flatten()[g_idx]

        for p_idx, policy in enumerate(POLICIES):

            policy_in_pub = PAPER_POLICY_LOOKUP[policy]

            data = experiment_data[policy, award_amount, '0.00', '0.00']

            fdr = data['falseDiscoveryRate'][:]
            fdr[np.isnan(fdr)] = 0.0

            fpr = data['falsePositiveRate'][:]

            # To 100 to plot only first 1M of 10M iterations.
            ax.plot(fdr.mean(axis=0)[:100],
                    label='{}, $F$'.format(policy_in_pub),
                     color=lcs[p_idx])
            ax.plot(fpr.mean(axis=0)[:100],
                    label=r'{}, $\overline{{\alpha}}$'.format(policy_in_pub),
                     color=lcs[p_idx], ls='--')

        ax.set_title('$G={}$'.format(award_amount))
        ax.set_xticks([0, 50, 100])
        ax.set_xticklabels(['0', '5e5', '1e6'])

        ax.set_yticks(np.arange(0, 1.01, 0.25))
        if g_idx % 2 == 0:
            ax.set_ylabel(r'$\overline{\alpha}$, $F$', size=12)
        if g_idx > 1:
            ax.set_xlabel('Iteration', size=12)
        if g_idx == 3:
            ax.legend(ncol=2)  #, bbox_to_anchor=(1.5, 0.4))
            ax.legend(handlelength=1, ncol=3, bbox_to_anchor=(.08, 0.4))

        ax.grid(True)

    if save_path is not None:
        plt.savefig(save_path)


def fpr_allpi(experiment_data, param='NPR', G='10', figsize=(8, 1.5),
              nSteps=100, nAgents=100, lcs=['blue', 'red', 'black'],
              param_vals=['0.10', '0.50', '0.90'],
              save_path=None):

    ax_idx = 0
    fig, axes = plt.subplots(nrows=1, ncols=3, figsize=figsize, sharey=True)

    x = np.array([jj + 1 for ii in range(nSteps // 5)
                     for jj in [ii]*nAgents
                 ])

    policies = ['PUBLICATIONS', 'RANDOM', 'FPR']
    policies_labels = ['PH', 'RA', 'MI']

    def _draw_legend(bbox_to_anchor=(0.15, 0.1)):

        lines = [
            mlines.Line2D([], [], color=lcs[policy_idx],
            markersize=10, marker='.', lw=0, label=policy)
            for policy_idx, policy in enumerate(policies_labels)
        ]
        lg = ax.legend(handles=lines, loc='lower left',
                       bbox_to_anchor=bbox_to_anchor)

    for pv_idx, pv in enumerate(param_vals):
        ax = axes[ax_idx]
        for policy_idx, policy in enumerate(policies):

            if param == 'NPR':
                try:
                    d = experiment_data[policy, G, pv, '0.00']
                except:
                    print(policy, G, pv)
                    print(list(experiment_data[policy, G, pv].keys()))
            else:
                d = experiment_data[policy, G, '0.00', pv]

            y = d['agentFPRs'][:].flatten()[::5]

            ax.plot(x, y, '.', color=lcs[policy_idx], ms=1.5, label=policy)
            if G != 60:
                if ax_idx == 0:
                    if G == 35 and param == 'FPDR' and pv == '0.10':
                        _draw_legend(bbox_to_anchor=(0.15, 0.2))
                    else:
                        _draw_legend()

            else:
                if param == 'NPR' and ax_idx == 1 and pv == '0.50':
                    _draw_legend()
                elif param == 'FPDR' and ax_idx == 2 and pv == '0.90':
                    _draw_legend(bbox_to_anchor=(0.25, 0.1))

        # if pv_idx == 0:
        if param == 'NPR':
            param_label = r'$p$'
        else:
            param_label = r'$r$'
        ax.set_title(r'{}={}'.format(param_label, pv))
        ax.set_ylabel(r'$\alpha_i$', size=16)
        # else:
        #     ax.set_title(str(pv))

        ax.set_xlabel('Iteration', size=14)
        ax.set_xticks([1, 10, 20])
        ax.set_xticklabels(['1e3', '5e6', '1e7'])
        ax.grid(True)
        ax.set_xlim(0, 21)

        ax_idx += 1

    if save_path is not None:
        plt.savefig(save_path)


def supplemental_policy_heatmaps(
        json_dict=None,
        experiment_data_dir=os.path.join('data', 'scimod-mixed-policy'),
        policy='MIXED',
        measure='falsePositiveRate',
        award_amounts=[10, 35, 60, 85],
        fpdr_npr_rates=[0.0, 0.25, 0.5, 0.75, 1.0],
        policyParams=[0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    ):
    '''
    Because the metadata for policy parameter is not in the HDF, we have to
    read directly from the data directory. This isn't too big a deal since
    there are only 120 JSONs for this experiment per policy. Two policies:
    MIXED and MODIFIED_RANDOM.

    Makes eight heatmaps for each policy. One for average false positive rate
    and one for false discovery rate.
    '''
    if json_dict is None:
        jsons = (
            json.load(open(g))
            for g in glob(os.path.join(experiment_data_dir, '*.json'))
        )
        json_dict = {}
        for j in jsons:
            params = j['metadata']['parameters']
            award_amount = params['awardAmount']
            fpdr_npr = params['falsePositiveDetectionRate']
            policy_param = params['policyParam']
            json_dict.update({(award_amount, fpdr_npr, policy_param): j})

    for amt in award_amounts:
        amts = {k[1:]: v for k, v in json_dict.items() if k[0] == amt}
        #### PUT HEATMAP MAKING HERE, YO! ####
        x = fpdr_npr_rates
        y = policyParams

        data = np.zeros((len(y), len(x)))
        for x_i, x_el in enumerate(x):
            for y_i, y_el in enumerate(y):
                try:
                    predata = np.array(amts[x_el, y_el][measure])[:, -1]
                    predata[predata == None] = 0.0
                    data[y_i, x_i] = predata.mean()
                except:
                    import ipdb
                    ipdb.set_trace()

        plt.figure()
        ax = sns.heatmap(data, vmin=0, vmax=1, xticklabels=x, yticklabels=y, cmap='viridis')
        ax.invert_yaxis()

        plt.savefig('policy_heatmap-G={}.pdf'.format(amt))

    return json_dict, ax
