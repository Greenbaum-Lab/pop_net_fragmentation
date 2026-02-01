"""
centrality.py

This module provides functions to compute node centralities (degree and betweenness)
for networks and collections of networks, especially in the context of network fragmentation
analyses. It uses NetworkX for graph processing and pandas for data handling. The main
functions allow centrality calculations for a single network, across replicates (e.g., 
simulation steps), and across multiple fragmentation types.
"""

from typing import List, Dict, Literal

import networkx as nx
import pandas as pd
from scipy.stats import pearsonr

from funcs import FragmentationResult


def compute_centrality_network(graph: nx.Graph) -> pd.DataFrame:
    """
    Compute degree and betweenness centrality for all nodes in a network.

    :param graph: NetworkX graph instance.
    :return: DataFrame with columns ['node_number', 'degree_centrality', 'betweenness_centrality'].
    """
    # degree_centrality = nx.degree_centrality(graph)
    betweenness_centrality = nx.betweenness_centrality(graph)
    degree_centrality = dict(nx.degree(graph))

    #         'degree': lambda net: dict(nx.degree(net))

    df = pd.DataFrame({
        'node_number': list(degree_centrality.keys()),
        'degree': list(degree_centrality.values()),
        'betweenness': list(betweenness_centrality.values())
    })

    return df


def compute_centrality_replicates(
    networks: List[List[nx.Graph]]
) -> pd.DataFrame:
    """
    Compute node centralities for all replicate-step networks.

    :param networks: Nested list of graphs [replicate][step].
    :return: DataFrame with columns ['replica', 'step', 'node_number', 'degree', 'betweenness'].
    """
    records = []

    for replica_idx, replicate_graphs in enumerate(networks):
        for step_idx, graph in enumerate(replicate_graphs):
            centralities_df = compute_centrality_network(graph)
            centralities_df['replica'] = replica_idx
            centralities_df['step'] = step_idx
            records.append(centralities_df)

    return pd.concat(records, ignore_index=True)


def compute_centrality_types(
    data: Dict[str, FragmentationResult],
    frag_types: list[str]
) -> pd.DataFrame:
    """
    Compute degree and betweenness centralities for all graphs across multiple fragmentation types.

    :param data: Mapping from frag_type to FragmentationResult.
    :param frag_types: List of fragmentation type keys to process.
    :return: DataFrame with columns ['fragmentation_type', 'replica', 'step', 'node_number', 'degree_centrality', 'betweenness_centrality'].
    """
    all_dfs = []
    for frag_type in frag_types:
        frag_res = data[frag_type]
        df = compute_centrality_replicates(frag_res.networks)
        df['frag_type'] = frag_type
        all_dfs.append(df)

    combined_df = pd.concat(all_dfs, ignore_index=True)
    cols = ['frag_type', 'replica', 'step', 'node_number', 'degree', 'betweenness']
    combined_df.to_csv(f'./csv_new/centrality.csv', index=False)
    return combined_df[cols]

