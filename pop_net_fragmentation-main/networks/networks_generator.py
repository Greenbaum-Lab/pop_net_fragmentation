from __future__ import annotations

import itertools as it
from typing import Iterable, List, Optional

import networkx as nx
import numpy as np

from Transformation import conservative_from_normal


# --------------------------------------------------------------------- #
# 1.  Undirected → directed (asymmetric) conversion
# --------------------------------------------------------------------- #

def symm_to_asymm_net(
    net: nx.Graph,
    mu: float = 1.0,
    sigma: float = 0.3,
    lower: float = 0.1,
    upper: float = 4,
    rng: Optional[np.random.Generator] = None,
) -> nx.DiGraph:
    """
    Replace every undirected edge in *net* by **two** directed edges whose
    weights are drawn from 𝒩(μ, σ²) and then mass-balanced so that

        Σⱼ mᵢⱼ  =  Σⱼ mⱼᵢ      for every node i.              (conservative)

    Parameters
    ----------
    net
        Undirected NetworkX graph carrying optional node / graph attrs.
    mu, sigma
        Parameters of the normal proposal before balancing.
    lower, upper
        Clip weights to `[lower, upper]` **after** balancing.
    rng
        Optional `numpy.random.Generator` for reproducibility.  If *None*
        we fall back to NumPy’s global PRNG.

    Returns
    -------
    nx.DiGraph
        New directed graph with identical nodes & graph attributes.
    """
    rng = np.random.default_rng() if rng is None else rng

    # (1) support matrix (binary, symmetric)
    support = nx.to_numpy_array(net, weight=None, dtype=int)

    # (2) conservative directed matrix on that support
    directed = conservative_from_normal(
        support,
        mu=mu,
        sigma=sigma,
        lower=lower,
        upper=upper,
        seed=int(rng.integers(2**32 - 1)),
    )

    # (3) build output DiGraph, preserving metadata
    out = nx.DiGraph()
    out.add_nodes_from((v, d) for v, d in net.nodes(data=True))
    out.graph.update(net.graph)

    # (4) copy edge weights (speed: vectorised lookup via dict)
    index = {node: i for i, node in enumerate(net.nodes())}
    for u, v in net.edges():
        i, j = index[u], index[v]
        out.add_edge(u, v, weight=float(directed[i, j]))
        out.add_edge(v, u, weight=float(directed[j, i]))

    return out


# --------------------------------------------------------------------- #
# 2.  Random-geometric graphs (RGG)
# --------------------------------------------------------------------- #

def make_rgg(
    n_graphs: int,
    n_nodes: int,
    target_edges: int,
    *,
    asymmetric: bool,
    radius: float = 0.30,
    rng: Optional[np.random.Generator] = None,
) -> List[nx.Graph | nx.DiGraph]:
    """
    Draw *n_graphs* random-geometric graphs until each has *target_edges*.

    Parameters
    ----------
    n_graphs, n_nodes, target_edges
        Self-explanatory.
    asymmetric
        If ``True`` convert each undirected RGG via :func:`symm_to_asymm_net`.
    radius
        Euclidean radius passed to :func:`nx.random_geometric_graph`.
    rng
        Optional `numpy.random.Generator` for reproducible draws.

    Returns
    -------
    list of graphs
    """
    rng = np.random.default_rng() if rng is None else rng
    graphs: List[nx.Graph | nx.DiGraph] = []

    while len(graphs) < n_graphs:
        g = nx.random_geometric_graph(n_nodes, radius, seed=int(rng.integers(2**32 - 1)))
        if g.number_of_edges() == target_edges:
            graphs.append(symm_to_asymm_net(g, rng=rng) if asymmetric else g)

    return graphs


# --------------------------------------------------------------------- #
# 3.  Spatial Erdős–Rényi (ER) helpers
# --------------------------------------------------------------------- #

def make_spatial_er(n_nodes: int, p: float, rng: Optional[np.random.Generator] = None) -> nx.Graph:
    """Single ER graph with random 2-D coordinates in node attr ``'pos'``."""
    rng = np.random.default_rng() if rng is None else rng
    positions = {i: rng.random(2) for i in range(n_nodes)}
    g = nx.erdos_renyi_graph(n_nodes, p, seed=int(rng.integers(2**32 - 1)))
    nx.set_node_attributes(g, positions, "pos")
    return g


def make_spatial_er_nets(
    n_graphs: int, n_nodes: int, p: float, rng: Optional[np.random.Generator] = None
) -> List[nx.Graph]:
    """Return *n_graphs* independent spatial ER networks."""
    rng = np.random.default_rng() if rng is None else rng
    return [make_spatial_er(n_nodes, p, rng) for _ in range(n_graphs)]


# --------------------------------------------------------------------- #
# 4.  Spatial small-world (SW) square-grid helpers
# --------------------------------------------------------------------- #

def make_spatial_sw(dim: int = 2, p: float = 0.015, rng: Optional[np.random.Generator] = None) -> nx.Graph:
    """
    Square grid (no periodic wrap) + random long-range rewiring.

    The grid has ``dim × dim`` nodes.  Each non-diagonal pair (i, j) is
    connected with probability *p*.
    """
    rng = np.random.default_rng() if rng is None else rng

    g = nx.grid_graph(dim=[dim, dim], periodic=False)
    mapping = {node: i for i, node in enumerate(g.nodes())}
    g = nx.relabel_nodes(g, mapping)               # compact labels 0..n-1

    # store coordinates for plotting later
    coords = {new: old for old, new in mapping.items()}
    nx.set_node_attributes(g, coords, "pos")

    nodes = list(g.nodes())
    for i, j in it.combinations(nodes, 2):
        if rng.random() < p:
            g.add_edge(i, j)

    return g


def make_spatial_sw_nets(
    n_graphs: int, dim: int = 2, p: float = 0.015, rng: Optional[np.random.Generator] = None
) -> List[nx.Graph]:
    """Return *n_graphs* independent spatial small-world networks."""
    rng = np.random.default_rng() if rng is None else rng
    return [make_spatial_sw(dim, p, rng) for _ in range(n_graphs)]


# --------------------------------------------------------------------- #
# 5.  Projection to conservative migration (re-balancing)
# --------------------------------------------------------------------- #


def project_to_conservative(
    net: nx.Graph,
    lower: float = 0.1,
    upper: float = 4,
    rtol: float = 1e-10,
) -> nx.Graph:
    """
    Adjust edge weights *in-place* so every weakly-connected component
    satisfies the conservative condition (row-sum == col-sum).

    Parameters
    ----------
    net
        Directed **or** undirected graph with ``weight`` attributes.
    lower, upper
        Optional clipping bounds for edge weights *after* balancing.
    rtol
        Relative tolerance passed to the Moore–Penrose pseudo-inverse.
    """
    # pick connected-component iterator depending on graph type
    comps: Iterable[set[int]] = (
        nx.weakly_connected_components(net) if net.is_directed()
        else nx.connected_components(net)
    )

    for comp in comps:
        if len(comp) <= 1:
            continue

        nodes = list(comp)
        M = nx.to_numpy_array(net, nodelist=nodes, weight="weight", nonedge=0.0)

        nz_idx = np.argwhere(M != 0)          # each row → edge index
        if nz_idx.size == 0:
            continue

        w = M[nz_idx[:, 0], nz_idx[:, 1]].astype(float)

        n, m = len(nodes), len(w)
        A = np.zeros((n, m))
        for k, (i, j) in enumerate(nz_idx):
            A[i, k] = 1.0
            A[j, k] = -1.0

        # ――― Quick-patch: compute pseudo-inverse *once* ――― #
        pinv = np.linalg.pinv(A @ A.T, rcond=rtol)

        # first balancing pass
        imbalance = A @ w
        if not np.allclose(imbalance, 0, atol=1e-12):
            w -= A.T @ pinv @ imbalance

        # optional clipping, then second (rare) balancing pass
        if lower is not None or upper is not None:
            w = np.clip(w, lower, upper)
            imbalance = A @ w
            if not np.allclose(imbalance, 0, atol=1e-12):
                w -= A.T @ pinv @ imbalance

        # write weights back to graph
        for val, (i, j) in zip(w, nz_idx):
            u, v = nodes[i], nodes[j]
            net[u][v]["weight"] = float(val)

    return net
