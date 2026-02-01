from __future__ import annotations

from typing import Iterable, Sequence, Tuple

import numpy as np
import pandas as pd
import networkx as nx


# --------------------------------------------------------------------- #
# Internal helpers
# --------------------------------------------------------------------- #

def _pairwise_flat(matrix: np.ndarray) -> np.ndarray:
    """
    Return the upper-triangle (i < j) values of *matrix* as a 1-D array.
    """
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("pairwise flatten expects a square 2-D matrix")

    inds = np.triu_indices_from(matrix, k=1)
    return matrix[inds]


def _assemble_distribution(values: Sequence[np.ndarray], genetics: str) -> pd.DataFrame:
    """Stack 1-D arrays into tidy long-format DataFrame."""
    records = [
        (step, v)
        for step, arr in enumerate(values)
        for v in _pairwise_flat(arr)
    ]
    return pd.DataFrame(records, columns=["step", genetics])


def _summarise(df: pd.DataFrame, genetics: str) -> pd.DataFrame:
    """Return mean of genetics per fragmentation step."""
    return (
        df.groupby("step")[genetics]
          .agg(avg="mean")
          .reset_index()
    )


# --------------------------------------------------------------------- #
# Public API
# --------------------------------------------------------------------- #

def calculate_genetics(
    migrations: Iterable[nx.Graph],
) -> Tuple[list[np.ndarray], list[np.ndarray]]:
    """
    Transform a list of migration graphs into per-step heterozygosity
    vectors and F_ST matrices.

    Returns
    -------
    het_list, fst_list
        • het_list – list[np.ndarray]  (shape = n_nodes) per step
        • fst_list – list[np.ndarray]  (shape = n_nodes × n_nodes) per step
    """
    from Transformation import transform_matrix

    het_list, fst_list = [], []

    for graph in migrations:
        M = nx.to_numpy_array(graph, weight="weight")
        T, F = transform_matrix(M)

        het = np.diag(T) / len(M)         # expected heterozygosity
        het_list.append(het.astype(float))
        fst_list.append(F.astype(float))

    return het_list, fst_list


def make_fst_dist(fst_mats: Sequence[np.ndarray]) -> pd.DataFrame:
    """Flatten a sequence of F_ST matrices into a tidy DataFrame."""
    return _assemble_distribution(fst_mats, genetics="fst")


def make_het_dist(het_vecs: Sequence[np.ndarray]) -> pd.DataFrame:
    """
    Flatten heterozygosity vectors into a tidy DataFrame.

    Columns
    -------
    • step  – fragmentation index
    • het   – heterozygosity of a single population
    """
    return (
        pd.DataFrame(het_vecs)
          .stack()                       # long format
          .rename_axis(("step", "pop"))
          .reset_index(name="het")
          .drop(columns="pop")
    )


def make_fst_stat(fst_df: pd.DataFrame) -> pd.DataFrame:
    """Mean F_ST per step."""
    return _summarise(fst_df, "fst")


def make_het_stat(het_df: pd.DataFrame) -> pd.DataFrame:
    """Mean heterozygosity per step."""
    return _summarise(het_df, "het")
