import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns

from distance_matrices import (get_random_walk_matrix, get_euclidean_matrix,
                               get_shortest_path_matrix, get_resistance, get_mfpt)
from funcs import load_data
from cor_test import calculate_qap


def plot_mantel(df):
    """
    plot mantel correlation for all fragmentation types.
    """
    df['step'] = df['step'] / df['step'].max() * 100  # Normalize steps to percentage
    df = df[df['p_val'] < 0.05]  # Filter rows with pval < 0.05
    # Filter steps with at least 5 unique replicas
    df = (
        df
        .groupby(['fragmentation_type', 'step'])
        .filter(lambda g: g['replica'].nunique() >= 5)
    )
    plt.figure(figsize=(10, 6))
    sns.lineplot(x='step', y='r_val', hue='fragmentation_type', data=df, errorbar='sd', legend=False)
    plt.xlabel('% fragmentation', fontsize=28)
    plt.ylabel('Correlation (r)', fontsize=28)
    plt.tick_params(axis='both', labelsize=25)
    plt.ylim(-0.05, 1.1)
    plt.savefig(f'fst_path.svg', format="svg")
    plt.show()



## plot single correlation fst-distance
# fragmentation_types = ['wrst']
# data = load_data(fragmentation_types)
# data = data[fragmentation_types[0]]
#
# steps = [0, 75, 150]
# fig, axes = plt.subplots(1, 3, figsize=(18, 7), sharey=True)
#
# for i, step in enumerate(steps):
#     net = data[1][5][step]
#     fst = data[7][5][step]
#     distance_matrix = get_random_walk_matrix(net)
#     r, p = calculate_mantel(net=net, fst_matrix=fst, dist_type='random', perms=999)
#     print(r, p)
#
#     flat_matrix1 = distance_matrix.flatten()
#     flat_matrix2 = fst.flatten()
#     flat_matrix1 = flat_matrix1[flat_matrix1 != 0]
#     flat_matrix2 = flat_matrix2[flat_matrix2 != 0]
#     df = pd.DataFrame({'distance': flat_matrix1, 'fst': flat_matrix2})
#     df = df.dropna()
#     # filter inf
#     df = df[~df['distance'].isin([np.inf, -np.inf])]
#     print(df)
#
#     sns.regplot(x='distance', y='fst', data=df, fit_reg=True, order=1, ax=axes[i])
#     axes[i].set_xlabel('Distance', fontsize=30)
#     axes[i].set_ylabel(r'Pairwise $F_{ST}$' if i == 0 else '', fontsize=30)
#     axes[i].tick_params(axis='both', labelsize=25)
#     axes[i].set_ylim(0, 0.5)
#     axes[i].text(0.05, 1.2, f'r={r:.2f}\np={p:.2e}', fontsize=20, transform=axes[i].transAxes)
#
# plt.tight_layout()
# plt.savefig(f'./figs/random_fst_steps.svg', format="svg")
# plt.show()



def plot_fst_distance(fragmentation_type, steps, dist_type, perms=999, output_path='test.svg'):
    """
    Plot single correlation FST-distance for specified fragmentation types and steps.
    :param fragmentation_type: List of fragmentation types to process
    :param steps: List of steps to process
    :param dist_type: Distance type identifier ('random', 'euclidean', 'path', 'resistance')
    :param perms: Number of permutations for Mantel test
    :param output_path: Path to save the output plot
    """
    data = load_data(fragmentation_type)
    data = data[fragmentation_type[0]]

    fig, axes = plt.subplots(1, len(steps), figsize=(18, 7), sharey=True)

    # choose distance matrix function
    dist_func_map = {
        'random': get_mfpt,
        'euclidean': get_euclidean_matrix,
        'path': get_shortest_path_matrix,
        'resistance': get_resistance,
        'mfpt': get_mfpt,
    }
    if dist_type not in dist_func_map:
        raise ValueError(f"Unknown dist_type: {dist_type}")
    dist_func = dist_func_map[dist_type]

    for i, step in enumerate(steps):
        net = data.networks[5][step]
        fst = data.fst_matrices[5][step]
        distance_matrix = dist_func(net)
        r, p = calculate_qap(net=net, fst_matrix=fst, dist_type=dist_type, perms=perms)

        flat_matrix1 = distance_matrix.flatten()
        flat_matrix2 = fst.flatten()
        flat_matrix1 = flat_matrix1[flat_matrix1 != 0]
        flat_matrix2 = flat_matrix2[flat_matrix2 != 0]
        df = pd.DataFrame({'distance': flat_matrix1, 'fst': flat_matrix2})
        df = df.dropna()
        # filter inf
        df = df[~df['distance'].isin([np.inf, -np.inf])]
        print(df)

        sns.regplot(x='distance', y='fst', data=df, fit_reg=True, order=1, ax=axes[i])
        axes[i].set_xlabel('Distance', fontsize=30)
        axes[i].set_ylabel(r'Pairwise $F_{ST}$' if i == 0 else '', fontsize=30)
        axes[i].tick_params(axis='both', labelsize=25)
        axes[i].set_ylim(0, 0.5)
        axes[i].text(0.05, 1.2, f'r={r:.2f}\np={p:.2e}', fontsize=20, transform=axes[i].transAxes)

    plt.tight_layout()
    plt.savefig(output_path, format="svg")
    plt.show()


