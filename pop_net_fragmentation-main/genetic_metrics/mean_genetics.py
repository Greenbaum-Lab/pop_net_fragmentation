"""
This module provides functionality for processing and visualizing genetic data
across different fragmentation types. It includes functions to normalize and
combine replicate-level mean data and plot mean ± SD of heterozygosity ('het')
or fixation index ('fst') across fragmentation types.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict

from funcs import FragmentationResult, percent_step


def mean_het_fst(data: Dict[str, FragmentationResult], measure: str) -> pd.DataFrame:
    """
    Combine and normalize replicate-level mean data for all fragmentation types.

    :param data: Mapping from frag_type to FragmentationResult.
    :param measure: 'het' or 'fst'.
    :return: DataFrame with columns ['step_pct', 'avg', 'replica', 'frag_type'].
    """
    if measure not in {'het', 'fst'}:
        raise ValueError(f"Invalid measure {measure!r}. Expected 'het' or 'fst'.")

    all_types = []
    for frag_type, frag_res in data.items():
        if measure == 'het':
            df = frag_res.het_mean.copy()
        else:  # measure == 'fst'
            df = frag_res.fst_mean.copy()

        df = percent_step(df, step_col='step', pct_col='step_pct')
        df['frag_type'] = frag_type
        all_types.append(df[['step_pct', 'avg', 'replica', 'frag_type']])

    return pd.concat(all_types, ignore_index=True)


def plot_genetics(data: Dict[str, FragmentationResult], measure: str) -> None:
    """
    Plot mean ± SD of the specified measure across all fragmentation types.

    :param data: Mapping from frag_type to FragmentationResult.
    :param measure: 'het' or 'fst'.
    """
    if measure not in {'het', 'fst'}:
        raise ValueError(f"Invalid measure {measure!r}. Expected 'het' or 'fst'.")

    # Process all fragmentation types to get a unified DataFrame
    df = mean_het_fst(data, measure)

    plt.figure(figsize=(10, 6))
    sns.lineplot(
        data=df,
        x='step_pct',
        y='avg',
        hue='frag_type',
        estimator='mean',
        ci='sd'
    )
    plt.xlabel('% fragmentation', fontsize=30)
    plt.ylabel(measure.capitalize(), fontsize=30)
    plt.tick_params(axis='both', labelsize=25)
    plt.legend(title='Type', fontsize=15)
    plt.tight_layout()
    plt.show()