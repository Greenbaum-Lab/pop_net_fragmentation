"""
This module provides utilities for analyzing and visualizing the heterozygosity 
of individual nodes in population fragmentation simulations. It includes functions 
to randomly select node indices per replica, extract corresponding heterozygosity 
data, and generate plots for selected nodes across simulation steps.

Functions:
    - select_random_nodes: Randomly selects node indices per replica.
    - extract_nodes: Extracts heterozygosity data for selected nodes.
    - plot_het_nodes: Plots heterozygosity over time for the selected nodes.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from typing import Dict
from funcs import assign_node_numbers, percent_step

def select_random_nodes(
    df: pd.DataFrame,
    per_replica: int = 1,
    nodes_per_step: int = 50
) -> Dict[int, np.ndarray]:
    """
    For each replica, choose `per_replica` random node indices.

    Parameters:
        df (pd.DataFrame): DataFrame containing the heterozygosity data.
        per_replica (int): Number of random nodes to select per replica.
        nodes_per_step (int): Number of nodes in each step.

    Returns:
        Dict[int, np.ndarray]: Dictionary with replica ids as keys and arrays of selected node indices as values.
    """
    df = assign_node_numbers(df, nodes_per_step)
    selections = {}
    for rep, sub in df.groupby("replica"):
        n_nodes = sub["node_number"].nunique()
        picks = np.random.choice(n_nodes, min(per_replica, n_nodes), replace=False)
        selections[int(rep)] = picks
    return selections

def extract_nodes(
    df: pd.DataFrame,
    selections: Dict[int, np.ndarray],
    nodes_per_step: int = 50
) -> pd.DataFrame:
    """
    Extract the heterozygosity data of selected nodes for each replica and step.

    Parameters:
        df (pd.DataFrame): DataFrame with 'node_number', 'step', 'replica', and 'het' values.
        selections (Dict[int, np.ndarray]): Dictionary with replicas as keys and selected node indices as values.
        nodes_per_step (int): Number of nodes per step.

    Returns:
        pd.DataFrame: DataFrame containing only selected nodes, with a unique 'id' per node-replica.
    """
    out = []
    for rep, nodes in selections.items():
        sub = df[df["replica"] == rep]
        for node in nodes:
            node_df = sub[sub["node_number"] == node].copy()
            node_df["id"] = f"n{node}_r{rep}"
            out.append(node_df)
    return pd.concat(out, ignore_index=True).drop(columns=['replica', 'node_number'])

def plot_het_nodes(
    df: pd.DataFrame,
    n_nodes: int = 10,
) -> None:
    """
    Plot the heterozygosity for selected nodes across steps.

    Parameters:
        df (pd.DataFrame): DataFrame with 'step', 'id', and 'het' values.
        n_nodes (int): Number of nodes to plot (randomly sampled from available nodes).
    """
    node_ids = df['id'].unique()
    selected_nodes = np.random.choice(node_ids, min(n_nodes, len(node_ids)), replace=False)
    df = percent_step(df, step_col='step', pct_col='step_pct')

    fig, ax = plt.subplots(figsize=(6, 4))
    for node_id in selected_nodes:
        node_data = df[df['id'] == node_id]
        ax.plot(node_data['step_pct'], node_data['het'], color='grey', alpha=0.5)

    ax.set_xlabel('Time', fontsize=16)
    ax.set_ylabel('Heterozygosity', fontsize=16)
    ax.tick_params(axis='both', which='major', labelsize=14)
    plt.show()
