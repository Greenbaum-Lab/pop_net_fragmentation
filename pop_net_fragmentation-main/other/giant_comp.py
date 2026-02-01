"""
Analysis utilities for evaluating the relationship between network fragmentation,
giant component size, and heterozygosity metrics across different fragmentation types.
Includes tools for merging, binning, and plotting these statistics.
"""

from typing import Dict, List, Optional
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import stats

from funcs import FragmentationResult, giant_component_over_steps
from typing import Dict, List


def het_component(
    data: Dict[str, FragmentationResult],
    frag_type: str
) -> pd.DataFrame:
    """
    For one fragmentation type, grab its precomputed mean het per replica-step
    and its giant-component fraction, merge them, and drop any zero-component rows.

    :param data: Mapping frag_type → FragmentationResult
    :param frag_type: Key of the fragmentation type to process
    :return: DataFrame with columns
             ['replica','step','avg_het','component'], filtered to component > 0.
    """
    frag_res = data[frag_type]
    # Mean heterozygosity per replica-step
    het_rep = (
        frag_res.het_mean
    )

    # 2. Giant-component fraction per replica-step
    comp_df = giant_component_over_steps(frag_res.networks)

    # 3. Merge and drop zero-size
    merged = pd.merge(het_rep, comp_df, on=['replica', 'step'], how='inner')
    return merged[merged['component'] > 0].reset_index(drop=True)

def het_component_types(
    data: Dict[str, FragmentationResult],
    frag_types: List[str] = None
) -> Dict[str, pd.DataFrame]:
    """
    For each fragmentation type, compute the merged heterozygosity vs. giant-component
    """
    all_types = {
        ft: het_component(data, ft)
        for ft in frag_types
    }
    return all_types


def bin_het_component(
    df: pd.DataFrame,
    n_bins: int = 20
) -> pd.DataFrame:
    """
    Bin component fractions into n_bins and compute mean±sd of avg_het in each bin.
    """
    binned = pd.cut(df['component'], bins=n_bins)
    stats = (
        df
        .groupby(binned)['avg']
        .agg(mean_het='mean', sd_het='std').reset_index()
    )
    stats['component_mid'] = stats['component'].apply(lambda interval: interval.mid)
    return stats[['component_mid','mean_het','sd_het']]



def plot_het_component(
    data: Dict[str, FragmentationResult],
    frag_types: List[str] = None,
    n_bins: int = 20,
    output: str = './figs/het_component.svg'
):
    """
    For each fragmentation type, prepare het vs. component data, bin it,
    and plot mean ± SD heterozygosity against fraction in the largest component.
    """
    fig, ax = plt.subplots(figsize=(10, 6))
    palette = plt.get_cmap('tab10')

    for i, ft in enumerate(frag_types):
        merged = het_component(data, ft)
        # 2. Bin component fractions and compute mean±SD het
        binned = bin_het_component(merged, n_bins=n_bins)

        color = palette(i)
        ax.scatter(
            binned['component_mid'],
            binned['mean_het'],
            label=ft,
            color=color
        )
        ax.errorbar(
            binned['component_mid'],
            binned['mean_het'],
            yerr=binned['sd_het'],
            fmt='o',
            color=color,
            alpha=0.7
        )

    ax.set_xlabel('Fraction of nodes in largest component', fontsize=16)
    ax.set_ylabel('Heterozygosity', fontsize=16)
    ax.tick_params(labelsize=12)
    plt.tight_layout()
    plt.savefig(output, format='svg')
    plt.show()

