import warnings
from mpl_toolkits.axes_grid1 import AxesGrid
warnings.simplefilter('ignore')

import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle

from glob import glob
from collections import Counter


def plot_simple(dataPath):

    j = json.load(open(dataPath, 'r'))

    df = pd.DataFrame(
        {
            k: np.array(v).mean(axis=1)
            for k, v in j.items()
            if k != 'metadata'
        }
    )

    metadata = j['metadata']
    df.plot()

    return df


def plot_ten_fpr(dataDir, title=None, figsize=(4, 6), save_path=None):

    for f in glob(os.path.join(dataDir, "*.json")):

        j = json.load(open(f, 'r'))

        df = pd.DataFrame(
            {
                k: np.array(v).mean(axis=1)
                for k, v in j.items()
                if k != 'metadata'
            }
        )

        df.falsePositiveRate.plot()

    if title is not None:
        plt.title(title)


POLICIES = ['RANDOM', 'PUBLICATIONS']


def plot_means(data_dir='data', policies=POLICIES):

    results_df = pd.DataFrame(columns=policies)

    for policy in policies:
        # data_dir = os.path.join(root_data_dir, policy)
        trials = []
        for f in glob(os.path.join(data_dir, "*.json")):
            j = json.load(open(f, 'r'))
            # Calculate mean across PIs for each trial.
            trials.append(np.mean(j["falsePositiveRate"], axis=1))

        import ipdb
        ipdb.set_trace()
        results_df[policy] = np.mean(trials, axis=0)

    results_df.plot()

    return results_df


def _get_jsons(experiment_dir):
    return [json.load(open(f))
            for f in glob(os.path.join(experiment_dir, '*.json'))]


def _get_amount(json_):
    return json_['metadata']['parameters']['awardAmount']


def pub_plot_award_experiment(experiment_dir='../fundingExperiment',
                              save_path=None, fundings=[], figsize=(7.5, 8.5),
                              n_iter=1e6):

    assert len(fundings) == 6
    # Subplots dimensions.
    nrows = 3
    ncols = 2

    jsons = _get_jsons(experiment_dir)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=figsize,
                             sharex=True, sharey=True)

    print(Counter([j['metadata']['policy'] for j in jsons]))
    for r_idx in range(nrows):

        from_ = (ncols * r_idx)
        to = (ncols * (r_idx + 1))
        plot_award_experiment_row(
            jsons=jsons,
            fundings=fundings[from_:to],
            axes=axes[r_idx]
        )

    axes[2, 1].legend(loc='lower right', fontsize=12)

    # Customize axes.
    lab_size = 15
    axes[2, 0].set_xlabel(r'$t$', size=lab_size)
    axes[2, 1].set_xlabel(r'$t$', size=lab_size)

    for ax in axes[:, 0]:
        ax.set_ylabel(r'$\alpha$', rotation=0, size=lab_size)
        ax.yaxis.set_label_coords(-0.25, 0.5)

    for ax in axes.flatten():
        ax.set_yticks(np.arange(0.0, 1.01, 0.25))
        ax.set_xticklabels(ax.get_xticklabels(), size=lab_size - 3)
        ax.grid(axis='both')
        ax.set_xlim(0.0, 1e6)
        ax.set_ylim(-0.05, 1.05)
        ax.set_xticks([0.0, 5e5, 1e6])
        ax.set_xticklabels(['0', r'$5 \times 10^5$', r'$10^6$'],
                           size=lab_size - 3)

    if save_path is not None:
        plt.savefig(save_path)


class DataDir:

    def __init__(self, data_dir):

        self.data_dir = data_dir
        self.jsons = _get_jsons(data_dir)
        self.metadatas = [j['metadata'] for j in self.jsons]

        p = 'parameters'
        keys = (
            (m['policy'], m[p]['awardAmount'],
             m[p]['publishNegativeResultRate'],
             m[p]['falsePositiveDetectionRate'])
            for m in self.metadatas
        )

        self.policies = [m['policy'] for m in self.metadatas]

        self.fpdrs = [m[p]['falsePositiveDetectionRate']
                      for m in self.metadatas]

        self.pubneg_rates = [m[p]['publishNegativeResultRate']
                             for m in self.metadatas]

        self.awardAmounts = [m[p]['awardAmount'] for m in self.metadatas]

        self.d = dict(zip(keys, self.jsons))

    def __getitem__(self, key):

        return self.d[key]


def heatmaps(data, award_amounts=[10, 40, 70],
             policies=['FPR', 'RANDOM', 'PUBLICATIONS'],
             fpdrs=[0.1, 0.25, 0.75, 0.9],
             negres_pub_rates=[0.1, 0.25, 0.75, 0.9]):

    if type(data) is str:
        d = DataDir(data).d
    elif type(data) is dict:
        d = data

    # fig, axes = plt.subplots(
    #     nrows=len(award_amounts), ncols=len(policies),
    #     figsize=(8, 11.5)
    # )
    fig = plt.figure(figsize=(7, 10))
    grid = AxesGrid(fig, 111, nrows_ncols=(len(award_amounts), len(policies)),
                    axes_pad=.65, cbar_mode='single', cbar_location='bottom',
                    cbar_pad=0.25, cbar_size='4%'
                    )
    # fig.subplots_adjust(top=0.1) #, top=1.0)
    # cax = fig.add_axes([0.2, 0.01, 0.5, 0.04])

    mats = []
    ax_idx = 0
    for award_idx, award_amount in enumerate(award_amounts):
        for policy_idx, policy in enumerate(policies):
            hmdict = {
                (k[2], k[3]): _ave_final_pol(v)
                for k, v in d.items()
                if k[0] == policy and k[1] == award_amount
            }

            hmdata = np.zeros((4, 4))
            for ii, x in enumerate(fpdrs):
                for jj, y in enumerate(negres_pub_rates):
                    hmdata[ii, jj] = hmdict[(x, y)]

            ax = grid[ax_idx]
            ax_idx += 1
            mats.append(ax.imshow(hmdata, origin='lower', vmin=0.0, vmax=1.0))

            # if award_idx == len(award_amounts) - 1:
            ax.set_xlabel('False positive discovery rate')
            if policy_idx == 0:
                ax.set_ylabel('Neg. Res. Pub. Rate')  # , size=14)
            ax.xaxis.set_ticks_position('bottom')

            ax.set_xticklabels([''] + [fpdrs[0], fpdrs[2]], rotation=0)
            ax.set_yticklabels([''] + negres_pub_rates)

            ax.set_title('G={}, {}'.format(award_amount, policy))

    # fig.colorbar(mat, ax=cax, orientation='horizontal')
    # fig.colorbar(mats[-2], ax=axes[0,-1], orientation='vertical')
    # fig.colorbar(mats[-2], ax=axes[1,-1], orientation='vertical')
    # fig.colorbar(mats[-2], ax=axes[2,-1], orientation='vertical')
    # fig.colorbar(mats[-2], ax=axes[3,-1], orientation='vertical')
    cbar = grid.cbar_axes[0].colorbar(mats[0])
    cbar.ax.set_xlabel('Average false positive rate')
    cbar.ax.set_xticks([0.0, 0.5, 1.0])
    # plt.tight_layout()
    plt.savefig('/Users/mt/Desktop/test-heatmap.pdf', bbox_inches='tight')

    return d


def _ave_final_pol(j):
    return np.array(j['falsePositiveRate'])[:, -1].mean()


def plot_award_experiment_row(
            experiment_dir='fundingExperiment', jsons=None, save_path=None,
            fundings=None, low_funding=50, high_funding=100,
            figsize=(8, 2), axes=None, policy=None, n_iter=1e6,
            measure='falsePositiveRate'
        ):

    if jsons is None:
        jsons = _get_jsons(experiment_dir)

    # Need to define how to check for desired funding amounts to plot.
    if fundings is None:
        def _funding_check(funding):
            return (low_funding <= funding and funding <= high_funding)
    else:
        def _funding_check(funding):
            return (funding in fundings)

    # Build key=G, value=FPR-timeseries dictionary.
    if measure == 'falsePositiveRate':
        fpr_dict = {
            (_get_amount(j), j['metadata']['policy']):
            np.mean(j[measure], axis=0)
            for j in jsons
            if _funding_check(_get_amount(j))
        }
    elif measure == 'falseDiscoveryRate':
        fpr_dict = {
            (_get_amount(j), j['metadata']['policy']):
            np.mean(j[measure], axis=0)
            for j in jsons
            if _funding_check(_get_amount(j))
        }

    # Ascending ordered list of amounts.
    amounts = list({_get_amount(j) for j in jsons})
    amounts.sort()
    amounts = [
        amount for amount in amounts
        if _funding_check(amount)
    ]
    n_amounts = len(amounts)

    # If axes have not been passed as arg, create them.
    if axes is None:
        _, axes = plt.subplots(nrows=1, ncols=n_amounts, figsize=figsize)

    # Ready plot parameters
    t = np.arange(0, n_iter, 2000)
    labels = ['Random', 'Publications', 'False pos. rate']
    colors = ['red', 'blue', 'black']

    # Plot FPR timeseries for different values of G, G increases from left
    # to right. Each subplot has three lines, one for each award policy.
    for idx, amount in enumerate(amounts):
        try:
            ax = axes[idx]
        except:
            print("excepted")
            ax = axes

        if policy is not None:
            ax.plot(t, fpr_dict[(amount, policy)], color='k',
                    label=policy)

        else:
            for p_idx, policy in enumerate(['RANDOM', 'PUBLICATIONS', 'FPR']):
                to_plot = fpr_dict[(amount, policy)]
                print(to_plot.shape)
                print((amount, policy))
                ax.plot(t, to_plot, color=colors[p_idx],
                        label=labels[p_idx])

        ax.set_title(r'$G={}$'.format(amount))

        ax.set_ylim(0.0, 1.0)
        ax.set_xticklabels(
                ['{:1.0e}'.format(el)
                 for idx, el in enumerate(ax.get_xticks())]
            )

    if save_path is not None:
        plt.savefig(save_path)


def alpha_v_G(experiment_dir='fundingExperiment', save_path=None,
              low_funding=35, high_funding=105, figsize=(4, 6)):

    sync_file = '.alpha_v_G_sync.pickle'

    if os.path.exists(sync_file):
        print(
            "Loading mean and stddev dictionaries from " + sync_file
        )
        with open(sync_file, 'rb') as sf:
            fpr_mean_dict, fpr_stddev_dict = pickle.load(sf)

    else:
        print(
            "Creating mean and stddev dictionaries, saving to " + sync_file
        )
        jsons = _get_jsons(experiment_dir)
        fpr_mean_dict = {
            _get_amount(j):
            np.mean(j['falsePositiveRate'], axis=0)[-1]
            for j in jsons
            if (j['metadata']['policy'] == 'FPR')
        }

        fpr_stddev_dict = {
            _get_amount(j):
            np.std(j['falsePositiveRate'], axis=0)[-1]
            for j in jsons
            if (j['metadata']['policy'] == 'FPR')
        }

        with open(sync_file, 'wb') as sf:
            pickle.dump([fpr_mean_dict, fpr_stddev_dict], sf)

    amounts = list(fpr_mean_dict.keys())
    amounts.sort()
    amounts = [
        el for el in amounts
        if low_funding <= el and el <= high_funding
    ]

    means = [fpr_mean_dict[amount] for amount in amounts]
    stddevs = [fpr_stddev_dict[amount] for amount in amounts]

    # plt.plot(amounts, means, 'o', color='black')

    plt.figure()

    plt.errorbar(amounts, means, yerr=stddevs, color='black', marker='o',
                 elinewidth=0.5)

    plt.yticks(np.arange(0, 1.1, 0.25))

    plt.xlabel(r'$G$', size=15)
    plt.ylabel(r'$\alpha$', size=15, rotation=0)
    plt.gca().yaxis.set_label_coords(-0.15, 0.5)
    plt.gca().grid(axis='both')

    if save_path:
        plt.savefig(save_path)

def _is_policy(j, policy):
    return policy == j['metadata']['policy']


def alpha_v_g_over_negrespub(
        experiment_dir='negativeResults/', save_path=None, low_funding=1,
        high_funding=105, figsize=(4, 6), negres_rates=[0.0, 0.1, 0.25, 0.5],
        policy='FPR'
        ):
    # sync_file = '.alpha_over_negrespub'
    plt.figure(figsize=figsize)

    # Need to make one plot for each neg result rate.
    jsons = _get_jsons(experiment_dir)

    def _is_negres_rate(j, negres_rate):
        return (negres_rate ==
                j['metadata']['parameters']['publishNegativeResultRate'])

    plot_colors = {
        'RANDOM': 'red',
        'PUBLICATIONS': 'blue',
        'FPR': 'black'
    }
    styles = ['-', ':', '-.', '--']
    for idx, negres_rate in enumerate(negres_rates):

        fpr_mean_dict = {
            _get_amount(j):
            np.mean(j['falsePositiveRate'], axis=0)[-1]
            for j in jsons
            if (_is_negres_rate(j, negres_rate) and _is_policy(j, policy))
        }

        fpr_stddev_dict = {
            _get_amount(j):
            np.std(j['falsePositiveRate'], axis=0)[-1]
            for j in jsons
            if (_is_negres_rate(j, negres_rate) and _is_policy(j, policy))
        }

        amounts = list(fpr_mean_dict.keys())
        amounts.sort()
        amounts = [
            el for el in amounts
            if low_funding <= el and el <= high_funding
        ]

        means = [fpr_mean_dict[amount] for amount in amounts]
        stddevs = [fpr_stddev_dict[amount] for amount in amounts]

        plt.errorbar(amounts, means, yerr=stddevs,
                     color=plot_colors[policy], marker='o',
                     ms=2, elinewidth=0.5, ls=styles[idx],
                     label='negresRate={:.2f}'.format(negres_rate))

    title = {
        'FPR': 'Least-FPR allocation policy',
        'RANDOM': 'Random allocation policy',
        'PUBLICATIONS': 'Most-publications allocation policy'
    }
    plt.title(title[policy])
    plt.xlabel('Award amount', size=14)
    plt.ylabel('Average final false positive rate', size=14)

    plt.yticks(np.arange(0, 1.01, 0.25))

    plt.gca().grid(axis='both')
    plt.gca().tick_params(axis='both', labelsize=12)

    plt.legend()

    if save_path is not None:
        plt.savefig(save_path)


def peer_review(
            experiment_dir='../data/peer-review', save_path=None,
            figsize=(10, 7.5), falsepos_detect_rates=[0.25, 0.5, 0.75, 1.0]
        ):
    fig, axes = plt.subplots(nrows=2, sharex=True, figsize=figsize)

    jsons = _get_jsons(experiment_dir)

    # false positive detection rate (fpdr)
    def _is_fpdr(j, fpdr):
        return (fpdr ==
                j['metadata']['parameters']['falsePositiveDetectionRate'])

    styles = ['-', ':', '-.', '--']
    colors = ['blue', 'red']
    titles = ['Most-publications allocation', 'Random allocation']
    for policy_idx, policy in enumerate(['PUBLICATIONS', 'RANDOM']):
        for fpdr_idx, fpdr in enumerate(falsepos_detect_rates):
            fpr_mean_dict = {
                _get_amount(j):
                np.mean(j['falsePositiveRate'], axis=0)[-1]
                for j in jsons
                if (_is_fpdr(j, fpdr) and _is_policy(j, policy))
            }

            fpr_stddev_dict = {
                _get_amount(j):
                np.std(j['falsePositiveRate'], axis=0)[-1]
                for j in jsons
                if (_is_fpdr(j, fpdr) and _is_policy(j, policy))
            }

            amounts = list(fpr_mean_dict.keys())
            amounts.sort()
            amounts = [
                el for el in amounts
                # if low_funding <= el and el <= high_funding
            ]

            means = [fpr_mean_dict[amount] for amount in amounts]
            stddevs = [fpr_stddev_dict[amount] for amount in amounts]

            axes[policy_idx].errorbar(
                amounts, means, yerr=stddevs, color=colors[policy_idx],
                marker='o', ms=2, elinewidth=0.5, ls=styles[fpdr_idx],
                label='fpdr={:.2f}'.format(fpdr)
            )

        axes[policy_idx].legend(fontsize=10)

        axes[policy_idx].set_ylabel('False positive rate', fontsize=14)
        axes[policy_idx].set_title(titles[policy_idx], fontsize=14)

        axes[policy_idx].set_yticks(np.arange(0.0, 1.01, 0.25))
        axes[policy_idx].grid(axis='both')

    axes[1].set_xlabel('Award amount', fontsize=14)
    if save_path is not None:
        plt.savefig(save_path)


def pubs_v_fpr(experiment_dir='../fundingExperiment',
               save_path=None, fundings=[10, 20, 50, 80, 90, 105],
               figsize=(6, 4)):
    jsons = _get_jsons(experiment_dir)

    # fig, axes = plt.subplots(ncols=len(fundings))

    def _funding_check(funding):
        return (funding in fundings)

    fpr_dict = {
        _get_amount(j): {
            'fpr': np.array(j['falsePositiveRate']),
            'pubs': np.array(j['nPublications'])
        }
        for j in jsons
        if _funding_check(_get_amount(j)) and j['metadata']['policy'] == 'FPR'
    }

    # Ascending ordered list of amounts.
    amounts = list({_get_amount(j) for j in jsons})
    amounts.sort()
    amounts = [
        amount for amount in amounts
        if _funding_check(amount)
    ]

    # for ax, amount in zip(axes, amounts):
    for amount in amounts:
        d = fpr_dict[amount]

        plt.plot(d['fpr'][:, -1], d['pubs'][:, -1], 'o',
                 label=r'$G={}$'.format(amount), alpha=0.5)

    plt.legend(fontsize=10, ncol=2)
    plt.xlabel(r'False positive rate at $t=T$', size=14)
    plt.ylabel(r'Publications at $t=T$', size=14)
    plt.title('Agency policy: least FPR', size=14)

    plt.gca().grid(axis='both')
    plt.gca().tick_params(axis='both', labelsize=12)

    if save_path is not None:
        plt.savefig(save_path)

    return fpr_dict
