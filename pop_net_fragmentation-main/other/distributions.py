# -----------------------------------------------------------------------------
# Utility functions for analyzing distributions of genetic diversity measures
# such as heterozygosity (het) and fixation index (fst) across different levels
# of network fragmentation.
# -----------------------------------------------------------------------------

import pandas as pd
import numpy as np
from typing import Literal, Tuple, List, Optional

from matplotlib import pyplot as plt
from funcs import FragmentationResult, percent_step


def filter_intervals(
    frag_res: FragmentationResult,
    measure: Literal['het', 'fst'],
    interval_pct: int = 25
) -> pd.DataFrame:
    """
    Select node-level measure data at fixed fragmentation-percent intervals.

    :param frag_res: One fragmentation result.
    :param measure: Column to filter ('het' or 'fst').
    :param interval_pct: Percentage spacing of intervals (must divide 100 evenly).
    :return: DataFrame with columns ['step_pct', 'replica', measure].
    """
    if measure not in {'het', 'fst'}:
        raise ValueError(f"Invalid measure {measure!r}. Expected 'het' or 'fst'.")

    # Pick the genetic data distribution
    df = frag_res.het_dist if measure == 'het' else frag_res.fst_dist

    # Compute continuous 0–100 step_pct
    df = percent_step(df, step_col='step', pct_col='step_pct')

    # Snap to nearest interval_pct multiple
    df['step_pct'] = (
        (df['step_pct'] / interval_pct)
        .round()
        .astype(int) * interval_pct
    )

    # Filter to only those snapped intervals
    allowed_intervals = set(range(0, 101, interval_pct))
    df_filtered = df[df['step_pct'].isin(allowed_intervals)].copy()

    return df_filtered[['step_pct', 'replica', measure]]


def compute_histogram(
    df: pd.DataFrame,
    measure: str,
    bins: int = 40
) -> Tuple[List[int], np.ndarray, List[np.ndarray]]:
    """
    Prepare histogram data for each step_pct layer.

    :param df: DataFrame with columns ['step_pct', measure].
    :param measure: Column to histogram ('het' or 'fst').
    :param bins: Number of bins for the histogram.
    :return: Tuple containing:
        - steps: Sorted unique step_pct values.
        - bin_edges: Array of bin edges.
        - hist_counts: List of count arrays for each step.
    """
    steps = sorted(df['step_pct'].unique(), reverse=True)
    hist_counts = []
    bin_edges = None

    for step in steps:
        values = df.loc[df['step_pct'] == step, measure].values
        counts, edges = np.histogram(values, bins=bins, density=True)
        hist_counts.append(counts)
        bin_edges = edges

    return steps, bin_edges, hist_counts


def plot_distribution(
    df: pd.DataFrame,
    measure: str,
    frag_type: str,
    bins: int = 40
) -> None:
    """
    Plot a ridgeline histogram of genetic diversity for one fragmentation type.

    :param df: DataFrame with columns ['step_pct', measure].
    :param measure: Column to plot ('het' or 'fst').
    :param frag_type: Identifier for the fragmentation type.
    :param bins: Number of bins for the histogram.
    """
    if measure not in {'het', 'fst'}:
        raise ValueError(f"Invalid measure {measure!r}. Expected 'het' or 'fst'.")

    # Compute histogram layers
    steps, bin_edges, hist_counts = compute_histogram(df, measure=measure, bins=bins)

    # Generate colormap based on measure
    n_steps = len(steps)
    cmap = (
        plt.get_cmap('YlGnBu')(np.linspace(0, 1, n_steps))
        if measure == 'het'
        else plt.get_cmap('YlOrRd')(np.linspace(0, 1, n_steps))
    )

    # Plot histogram layers
    fig, ax = plt.subplots(figsize=(4, 2 + 0.5 * n_steps))
    bin_width = bin_edges[1] - bin_edges[0]

    for i, (step, counts) in enumerate(zip(steps, hist_counts)):
        base = i * 6
        ax.bar(
            bin_edges[:-1],
            counts,
            width=bin_width,
            bottom=base,
            color=cmap[i],
            edgecolor='black',
            alpha=0.6,
            align='edge'
        )
        ax.hlines(base, bin_edges[0], bin_edges[-1], color='black', linewidth=0.5)

    # Customize plot appearance
    ax.set_yticks([])
    ax.set_xlabel(measure.capitalize(), fontsize=14)
    ax.set_xlim(bin_edges[0], bin_edges[-1])
    ax.set_ylim(0, 6 * n_steps + max(cnt.max() for cnt in hist_counts))
    ax.tick_params(axis='both', labelsize=12)
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)

    plt.title(f"{frag_type.capitalize()} Fragmentation", fontsize=16)
    plt.tight_layout()
    plt.show()