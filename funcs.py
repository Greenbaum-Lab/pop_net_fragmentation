from dataclasses import dataclass
from typing import Dict, List
import pickle
import os
import networkx as nx
import pandas as pd
import numpy as np
import logging


logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class FragmentationResult:
    """
    Container for results of a fragmentation simulation.
    """
    n_steps: int
    networks: List[nx.Graph]
    het_dist: pd.DataFrame
    het_mean: pd.DataFrame
    fst_dist: pd.DataFrame
    fst_mean: pd.DataFrame
    coalescence_list: List[np.ndarray]
    fst_matrices: List[np.ndarray]


def load_data(fragmentation_types: List[str]) -> Dict[str, FragmentationResult]:
    """
    Load fragmentation results for specified types, wrapping each in a dataclass.

    :param fragmentation_types: List of fragmentation type identifiers (e.g., ["cor", "rand"]).
    :param base_path: Directory where the files are stored.
    :param extension: File extension of the pickled results.
    :return: Dict mapping frag_type -> FragmentationResult.
    """
    results: Dict[str, FragmentationResult] = {}
    for ft in fragmentation_types:
        file_path =  f"RGG, {ft}_asymmetric.pickle"
        with open(file_path, "rb") as f:
            raw = pickle.load(f)
        results[ft] = FragmentationResult(
            n_steps         = raw[0],
            networks        = raw[1],
            het_dist        = raw[2],
            het_mean        = raw[3],
            fst_dist        = raw[4],
            fst_mean        = raw[5],
            coalescence_list= raw[6],
            fst_matrices    = raw[7],
        )
    logger.info(f"Finished loading fragmentation types: {fragmentation_types}")

    return results


def percent_step(df, step_col='step', pct_col='step_pct'):
    df[pct_col] = df[step_col] / df[step_col].max() * 100
    return df


def fraction_largest_component(
    G: nx.Graph
) -> float:
    """
    Compute the fraction of nodes in the largest connected component of G.

    :param G: A NetworkX graph.
    :return: Size of the largest connected component divided by total number of nodes.
    """
    # Find all components, pick the largest by cardinality
    largest = max(nx.connected_components(G), key=len)
    return len(largest) / G.number_of_nodes()


def giant_component_over_steps(
    all_nets: List[List[nx.Graph]]  # list of replicates, each a list of step‐graphs
) -> pd.DataFrame:
    """
    For each replicate and step, compute the fraction of nodes in the
    largest connected component.

    :param all_nets: nested list of networks: all_nets[replica][step] → Graph
    :return: DataFrame with columns ['replica', 'step', 'component']
    """
    records = []
    for replica, nets in enumerate(all_nets):
        for step, net in enumerate(nets):
            records.append({
                'replica': replica,
                'step':    step,
                'component': fraction_largest_component(net)
            })
    return pd.DataFrame.from_records(records)


def assign_node_numbers(df: pd.DataFrame, nodes_per_step: int = 50) -> pd.DataFrame:
    """
    Assigns node numbers for each node in each replica.

    :param df: DataFrame containing the heterozygosity data.
    :param nodes_per_step: Number of nodes in the network.
    :return: DataFrame with the original data and an additional 'node_number' column.
    """
    # Ensure that 'step' and 'replica' are in the DataFrame
    # Create an array of node numbers for each replica
    df['node_number'] = df.groupby('replica').cumcount() % nodes_per_step

    return df