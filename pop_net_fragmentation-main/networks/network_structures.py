import networkx as nx
import pandas as pd
from typing import List, Dict
from funcs import fraction_largest_component  # you already refactored this

def fraction_isolated_nodes(G: nx.Graph) -> float:
    """Fraction of nodes that are isolated."""
    return len(list(nx.isolates(G))) / G.number_of_nodes()

def fraction_secondary_components(G: nx.Graph, min_size: int = 4) -> float:
    """
    Fraction of nodes in components of size ≥ min_size,
    excluding the giant component.
    """
    giant = max(nx.connected_components(G), key=len)
    comps = [
        c for c in nx.connected_components(G)
        if c != giant and len(c) >= min_size
    ]
    return sum(len(c) for c in comps) / G.number_of_nodes()

def fraction_waste(G: nx.Graph, min_size: int = 2, max_size: int = 3) -> float:
    """
    Fraction of nodes in “waste” components: size between min_size and max_size.
    """
    comps = [c for c in nx.connected_components(G) if min_size <= len(c) <= max_size]
    return sum(len(c) for c in comps) / G.number_of_nodes()

def metrics_per_graph(G: nx.Graph) -> Dict[str, float]:
    """
    Compute all four fractions for a single graph.
    """
    return {
        'giant': fraction_largest_component(G),
        'isolated': fraction_isolated_nodes(G),
        'components': fraction_secondary_components(G),
        'waste': fraction_waste(G),
    }

def network_metrics_over_steps(
    networks: List[nx.Graph]
) -> pd.DataFrame:
    """
    For a single replicate (a list of graphs over steps), return a DataFrame:
    columns=['step','giant','isolated','components','waste'].
    """
    records = []
    for step, G in enumerate(networks):
        m = metrics_per_graph(G)
        m['step'] = step
        records.append(m)
    return pd.DataFrame.from_records(records)

def network_metrics_over_replicates(
    replicas: List[List[nx.Graph]]
) -> pd.DataFrame:
    """
    For multiple replicates, return a long DataFrame with columns
    ['replica','step','giant','isolated','components','waste'].
    """
    dfs = []
    for ridx, nets in enumerate(replicas):
        df_steps = network_metrics_over_steps(nets)
        df_steps['replica'] = ridx
        dfs.append(df_steps)
    return pd.concat(dfs, ignore_index=True)

def aggregate_network_metrics(
    df: pd.DataFrame
) -> pd.DataFrame:
    """
    Collapse over replicas by taking the mean of each metric at each step.
    Returns ['step','giant','isolated','components','waste'].
    """
    agg = (
        df
        .groupby('step')[['giant','isolated','components','waste']]
        .mean()
        .reset_index()
    )
    # ensure sums to 1 by recomputing 'waste' if you like:
    agg['waste'] = 1 - agg[['giant','isolated','components']].sum(axis=1)
    return agg



import matplotlib.pyplot as plt
from typing import Dict, List
from funcs import FragmentationResult




def plot_network_stacked_area_all(
    data: Dict[str, FragmentationResult],
    frag_types: List[str] = None,
    output_path: str = './figs/stack_all_fragmentation_types.svg'
):
    if frag_types is None:
        frag_types = list(data.keys())

    n = len(frag_types)
    cols = 2
    rows = (n + cols - 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(15, rows * 5), sharex=True, sharey=True)
    axes = axes.flatten()

    for i, ft in enumerate(frag_types):
        ax = axes[i]
        row_idx = i // cols
        col_idx = i % cols

        frag_res = data[ft]
        df_rep = network_metrics_over_replicates(frag_res.networks)
        df_agg = aggregate_network_metrics(df_rep)

        steps = df_agg['step']
        y = [df_agg[col] for col in ['waste','isolated','components','giant']]
        colors = plt.cm.Dark2.colors[:4]

        ax.stackplot(steps, *y, labels=['waste','isolated','components','giant'],
                     colors=colors, alpha=0.8)
        ax.set_title(ft, fontsize=18)
        ax.set_ylim(0,1)
        ax.tick_params(labelsize=14)

        # only bottom row gets x‐label
        if row_idx == rows - 1:
            ax.set_xlabel('Fragmentation step', fontsize=16)
        # only first column gets y‐label
        if col_idx == 0:
            ax.set_ylabel('Proportion', fontsize=16)

    # remove any unused subplots
    for j in range(n, len(axes)):
        fig.delaxes(axes[j])

    # common legend
    fig.legend(['waste','isolated','components','giant'],
               loc='upper center', ncol=4, frameon=False, fontsize=14)

    plt.tight_layout(rect=[0,0,1,0.95])
    plt.savefig(output_path, format='svg')
    plt.show()

from funcs import load_data


if __name__ == '__main__':
    # 1. load only the fragmentation types you care about
    frag_types = ['cor', 'rand', 'dist']
    data = load_data(frag_types)  # returns {ft: FragmentationResult}

    # 2. plot
    plot_network_stacked_area_all(
        data,
        frag_types=frag_types,
        output_path='./figs/stack_all_fragmentation_types.svg'
    )

