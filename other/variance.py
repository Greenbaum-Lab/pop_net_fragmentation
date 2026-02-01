"""
variance.py

Script to compute and plot the variance of heterozygosity for different fragmentation types
from simulation results. Functions robustly process a dictionary of results and return a
concatenated DataFrame with variance per replica and step for each specified fragmentation type.

"""

from typing import Dict, List
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os

from funcs import FragmentationResult, percent_step

def variance_per_replica_step(
    data: Dict[str, FragmentationResult],
    frag_type: str
) -> pd.DataFrame:
    """
    Calculate heterozygosity variance per (replica, step) for a given frag_type.

    :param data: Dict of fragmentation results.
    :param frag_type: Fragmentation type key.
    :return: DataFrame with columns ['replica', 'step', 'variance'].
    """
    frag_res = data[frag_type]
    df = frag_res.het_dist
    return (
        df.groupby(['replica', 'step'])['het']
          .var(ddof=1)
          .reset_index(name='variance')
    )


def process_variance(
    data: Dict[str, FragmentationResult],
    fragmentation_types: List[str]
) -> pd.DataFrame:
    """
    Prepare concatenated per-replica variance data.

    :param data: Dict of frag_type → FragmentationResult.
    :param fragmentation_types: List of frag_types to process.
    :return: DataFrame with columns ['fragmentation_type', 'replica', 'step', 'variance'].
    """
    dfs = []
    for frag_type in fragmentation_types:
        var_df = variance_per_replica_step(data, frag_type)
        var_df['frag_type'] = frag_type
        dfs.append(var_df)
    return pd.concat(dfs, ignore_index=True)

def plot_variance(df: pd.DataFrame) -> None:
    """
    Plot variance.

    :param df: DataFrame with columns ['fragmentation_type', 'replica', 'step', 'variance'].
    """
    plt.figure(figsize=(10, 6))
    # Normalize step to percentage
    df = percent_step(df, step_col='step', pct_col='step_pct')
    sns.lineplot(
        data=df,
        x='step_pct',
        y='variance',
        hue='frag_type',
        estimator='mean',
        ci='sd'
    )
    plt.xlabel('% fragmentation', fontsize=25)
    plt.ylabel('Variance', fontsize=25)
    plt.tick_params(axis='both', labelsize=20)
    plt.savefig('./figs/variance.svg', format='svg', dpi=300)
    plt.show()
