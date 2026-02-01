
"""
centrality_corr.py

This module provides functions for analyzing the correlation between network centrality measures 
and heterozygosity in network fragmentation analyses. Functions include merging centrality and 
heterozygosity data, computing correlations, filtering significant results, and plotting correlations.

"""
from typing import Dict
import pandas as pd
from funcs import FragmentationResult, assign_node_numbers
from matplotlib import pyplot as plt
from funcs import percent_step
from typing import Dict, Literal
from scipy.stats import pearsonr

def compute_het_central_correlation(
    df: pd.DataFrame,
    centrality: Literal['degree', 'betweenness'],
) -> pd.DataFrame:
    """
    Compute Pearson correlation (r) and p-value between centrality and heterozygosity
    for each (fragmentation_type, replica, step) group.

    :param df: DataFrame containing columns:
               ['frag_type', 'replica', 'step', centrality_col, heterozygosity_col]
    :param centrality_col: Name of centrality measure column ('degree' or 'betweenness')
    :param heterozygosity_col: Name of heterozygosity column (default 'het')
    :return: DataFrame with columns:
             ['frag_type', 'replica', 'step', 'r', 'p']
    """
    results = []

    # Ensure frag_type maintains its order
    frag_type_order = df['frag_type'].unique()
    df['frag_type'] = pd.Categorical(df['frag_type'], categories=frag_type_order, ordered=True)

    grouped = df.groupby(['frag_type', 'replica', 'step'])

    for (frag_type, replica, step), group in grouped:
        group = group[group[centrality] != 0]  # Exclude rows where centrality is 0
        if group[centrality].nunique() < 2:
            continue
        r, p = pearsonr(group[centrality], group['het'])
        results.append({
            'frag_type': frag_type,
            'replica': replica,
            'step': step,
            'r': r,
            'p': p
        })

    corr_df = pd.DataFrame(results)
    corr_df.to_csv(f'./csv_new/het_bet_correlation.csv', index=False)
    return pd.DataFrame(results)



def merge_centrality_het(
        centrality_df: pd.DataFrame,
        data: Dict[str, FragmentationResult],
        frag_types: list[str]
) -> pd.DataFrame:
    """
    Preprocess and merge the centrality data with heterozygosity data for each fragmentation type.

    :param centrality_df: DataFrame containing 'frag_type', 'replica', 'step', 'node_number', 'degree', 'betweenness'
    :param data: Dictionary mapping frag_type → FragmentationResult
    :param frag_types: List of fragmentation types to process
    :return: Merged DataFrame with centrality and heterozygosity for each node.
    """
    all_data = []

    # Iterate over each fragmentation type
    for frag_type in frag_types:
        # Get the heterozygosity data from FragmentationResult
        frag_res = data[frag_type]
        assign_node_numbers(frag_res.het_dist)
        het_df = frag_res.het_dist

        # Merge the centrality and heterozygosity data on ['replica', 'step', 'node_number']
        merged_df = pd.merge(
            centrality_df[centrality_df['frag_type'] == frag_type],
            het_df[['replica', 'step', 'node_number', 'het']],
            on=['replica', 'step', 'node_number'],
            how='left'  # 'left' join keeps all centrality data and adds 'het' where possible
        )

        all_data.append(merged_df)

    # Concatenate all fragmentation types into a single DataFrame
    final_df = pd.concat(all_data, ignore_index=True)
    final_df.to_csv(f'./csv_new/centrality_het.csv', index=False)

    return final_df



def filter_correlations(
    corr_df: pd.DataFrame,
    min_replicates: int
) -> pd.DataFrame:
    """
    Filter correlation DataFrame to include only significant results (p < threshold)
    and groups with more than min_replicates.

    :param corr_df: DataFrame with correlation results, including p-values.
    :param min_replicates: Minimum number of replicates required per (frag_type, step).
    :return: Filtered DataFrame.
    """
    df_filtered = corr_df[(corr_df['p'] < 0.05) & (corr_df['p'] > 0)]
    # Identify valid (frag_type, step) groups with enough replicates
    valid_groups = (
        df_filtered
        .groupby(['frag_type', 'step'])['replica']
        .nunique()
        .reset_index()
        .query(f"replica >= {min_replicates}")
        [['frag_type', 'step']]
    )

    return df_filtered.merge(valid_groups, on=['frag_type', 'step'], how='inner')




def plot_correlation(
    corr_df,
    output_path
, sns=None):
    """
    Plot correlation coefficient r over steps using Seaborn to compute mean ± SD.

    :param corr_df: DataFrame with columns ['frag_type', 'replica', 'step', 'r', 'p'].
    :param frag_type_col: Column name for fragmentation type.
    :param step_col: Column name for step.
    :param r_col: Column name for correlation coefficient.
    :param output_path: Path to save plot.
    """
    # Convert step to percentage using func percent_step
    corr_df = percent_step(corr_df, step_col='step', pct_col='step_pct')

    plt.figure(figsize=(6, 4))
    sns.lineplot(
        data=corr_df,
        x='step_pct',
        y='r',
        hue='frag_type',
        estimator='mean',
        errorbar='sd',
    )
    plt.xlabel('% fragmentation', fontsize=16)
    plt.ylabel('Correlation (r)', fontsize=16)
    plt.tick_params(axis='both', labelsize=14)
    plt.ylim(-1, 1.1)
    plt.legend().set_visible(False)
    plt.savefig(output_path, format='svg')
    plt.show()
