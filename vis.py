import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

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


def plot_mutation_parameters():
    pass
