from __future__ import annotations

import logging
from multiprocessing import Pool, cpu_count
from statistics import mean
from typing import Callable, List, Tuple

import networkx as nx
import numpy as np
import pandas as pd

from processes import (
    remove_edge_random,
    remove_edge_correlated,
    remove_edge_distance,
    remove_edge_intrusive,
    remove_edge_divisive,
    remove_edge_regressive,
    remove_edge_optimal,
    remove_edge_worst,
)

from stats_utils import (
    calculate_genetics,
    make_fst_dist,
    make_het_dist,
    make_fst_stat,
    make_het_stat,
)

from networks_generator import project_to_conservative  # keeps in-place contract

# ---------------------------------------------------------------------------- #
# 0.  Logging defaults
# ---------------------------------------------------------------------------- #

LOG = logging.getLogger(__name__)
if not LOG.handlers:
    h = logging.StreamHandler()
    h.setFormatter(logging.Formatter("%(levelname)s | %(message)s"))
    LOG.addHandler(h)
LOG.setLevel(logging.INFO)

# ---------------------------------------------------------------------------- #
# 1.  Mapping short tags → fragmentation functions
# ---------------------------------------------------------------------------- #

_FRAG_MAP: dict[str, Callable[[nx.Graph], List[nx.Graph]]] = {
    "rand": remove_edge_random,
    "cor":  remove_edge_correlated,
    "intr": remove_edge_intrusive,
    "reg":  remove_edge_regressive,
    "div":  remove_edge_divisive,
    "dist": remove_edge_distance,
    "opt":  remove_edge_optimal,
    "wrst": remove_edge_worst,
}

# ---------------------------------------------------------------------------- #
# 2.  Single-replicate runner
# ---------------------------------------------------------------------------- #

def run_single_fragmentation(
    net: nx.Graph,
    frag_key: str,
    replica: int,
) -> Tuple[
    int,                        # nets_number (snapshots per replicate)
    List[nx.Graph],             # list of snapshot graphs
    pd.DataFrame, pd.DataFrame, # het_dens, het_stat
    pd.DataFrame, pd.DataFrame, # fst_dens, fst_stat
    List[np.ndarray],           # genetics_coal
    List[np.ndarray],           # genetics_fst
]:
    """
    Apply one fragmentation scenario to *net* and compute genetic stats.
    """
    frag_fn = _FRAG_MAP[frag_key]

    # 1 ─ generate the migration sequence
    migration = frag_fn(G=net)
    migration = [project_to_conservative(g) for g in migration]

    n_snapshots = len(migration)

    # 2 ─ heavy algebra (coalescence + F_ST) – local, no I/O
    genetics_coal, genetics_fst = calculate_genetics(migration)

    # 3 ─ heterozygosity
    het_dens = make_het_dist(genetics_coal)
    het_dens["replica"] = replica

    het_stat = make_het_stat(het_dens)
    het_stat["replica"] = replica

    # 4 ─ F_ST
    fst_dens = make_fst_dist(genetics_fst)
    fst_dens["replica"] = replica

    fst_stat = make_fst_stat(fst_dens)
    fst_stat["replica"] = replica

    return (
        n_snapshots,
        migration,
        het_dens,
        het_stat,
        fst_dens,
        fst_stat,
        genetics_coal,
        genetics_fst,
    )

# ---------------------------------------------------------------------------- #
# 3.  Parallel replicates
# ---------------------------------------------------------------------------- #

def run_replicates(
    nets: List[nx.Graph],
    frag_key: str,
    n_workers: int | None = None,
) -> Tuple[
    float,                      # mean nets_number across replicates
    List[List[nx.Graph]],       # all_nets  (list-per-replicate)
    pd.DataFrame, pd.DataFrame, # het_dens, het_stat
    pd.DataFrame, pd.DataFrame, # fst_dens, fst_stat
    List[List[np.ndarray]],     # genetics_coal  (per replicate)
    List[List[np.ndarray]],     # genetics_fst
]:
    """
    Run `run_single_fragmentation` over all *nets* in parallel.

    Parameters
    ----------
    nets
        List of *base* migration graphs, each defining one replicate.
    frag_key
        Short tag in {'rand', 'cor', 'intr', …} – see `_FRAG_MAP`.
    n_workers
        Processes to spawn.  Default = `min(len(nets), cpu_count())`.

    Returns
    -------
    Tuple compatible with old `make_replicates_new`.
    """
    if frag_key not in _FRAG_MAP:
        raise KeyError(
            f"Unknown fragmentation key: {frag_key!r}. "
            f"Valid keys: {', '.join(_FRAG_MAP)}."
        )

    n_workers = min(len(nets), cpu_count()) if n_workers is None else n_workers
    LOG.info("Running %d replicates on %d worker(s)…", len(nets), n_workers)

    # --- fan out work ---------------------------------------------------- #
    with Pool(processes=n_workers) as pool:
        results = pool.starmap(
            run_single_fragmentation,
            [(net, frag_key, idx) for idx, net in enumerate(nets)],
        )

    # --- unpack results -------------------------------------------------- #
    (
        nets_number,    # tuple[float]
        all_nets,
        het_dens,
        het_stat,
        fst_dens,
        fst_stat,
        genetics_coal,
        genetics_fst,
    ) = zip(*results)

    # aggregate DataFrames
    het_dens = pd.concat(het_dens, ignore_index=True)
    het_stat = pd.concat(het_stat, ignore_index=True)
    fst_dens = pd.concat(fst_dens, ignore_index=True)
    fst_stat = pd.concat(fst_stat, ignore_index=True)

    return (
        mean(nets_number),
        list(all_nets),
        het_dens,
        het_stat,
        fst_dens,
        fst_stat,
        list(genetics_coal),
        list(genetics_fst),
    )

# ---------------------------------------------------------------------------- #
# 4.  Guard for Windows – allow “python -m pipeline” execution
# ---------------------------------------------------------------------------- #

if __name__ == "__main__":  # pragma: no cover
    import argparse
    from networks_generator import make_rgg  # simple demo

    parser = argparse.ArgumentParser(
        description="Quick smoke test of the fragmentation pipeline."
    )
    parser.add_argument("--replicates", type=int, default=4)
    parser.add_argument("--nodes", type=int, default=30)
    parser.add_argument("--edges", type=int, default=40)
    parser.add_argument("--frag", choices=_FRAG_MAP, default="rand")
    args = parser.parse_args()

    base_nets = make_rgg(
        n_graphs=args.replicates,
        n_nodes=args.nodes,
        target_edges=args.edges,
        asymmetric=False,
    )

    out = run_replicates(base_nets, frag_key=args.frag)
    LOG.info("Finished.  Heterozygosity stats head:\n%s", out[3].head())
