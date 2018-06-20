import warnings
warnings.simplefilter('ignore')

import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pickle

from glob import glob


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
    print(metadata)

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


def plot_award_experiment(experiment_dir='awardAmounts', save_path=None,
                          low_funding=50, high_funding=100, figsize=(8, 2)):

    jsons = _get_jsons(experiment_dir)

    fpr_dict = {
        (_get_amount(j), j['metadata']['policy']):
        np.mean(j['falsePositiveRate'], axis=0)
        for j in jsons
        if low_funding <= _get_amount(j) and _get_amount(j) <= high_funding
    }

    amounts = list({_get_amount(j) for j in jsons})
    amounts.sort()
    amounts = [
        el for el in amounts
        if low_funding <= el and el <= high_funding
    ]

    n_amounts = len(amounts)

    fig, axes = plt.subplots(nrows=1, ncols=n_amounts, figsize=figsize)

    t = np.arange(0, 1e6, 2000)
    labels = ['Random', 'Publications', 'False pos. rate']
    colors = ['red', 'blue', 'black']
    for idx, amount in enumerate(amounts):
        ax = axes[idx]
        for p_idx, policy in enumerate(['RANDOM', 'PUBLICATIONS', 'FPR']):
            ax.plot(t, fpr_dict[(amount, policy)], color=colors[p_idx],
                    label=labels[p_idx])

        ax.set_title(r'$G={}$'.format(amount))

        ax.set_ylim(0.0, 1.0)
        ax.set_xticklabels(
                ['{:1.0e}'.format(el)
                 for idx, el in enumerate(ax.get_xticks())]
            )

    if save_path is not None:
        plt.savefig(save_path)

    return jsons, fpr_dict


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

    plt.errorbar(amounts, means, yerr=stddevs, color='black', marker='o')

    plt.yticks(np.arange(0, 1.1, 0.25))

    plt.xlabel(r'$G$', size=15)
    plt.ylabel(r'$\alpha^*$', size=15, rotation=0)
    plt.gca().yaxis.set_label_coords(-0.15, 0.5)
    plt.gca().grid(axis='both')

    if save_path:
        plt.savefig(save_path)
