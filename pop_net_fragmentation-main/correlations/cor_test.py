import numpy as np
import networkx as nx
from typing import Optional, Tuple, List
import concurrent.futures
import pandas as pd
import pingouin as pg
from distance_matrices import get_mfpt
from scipy.stats import chi2
from distance_matrices import get_euclidean_matrix, get_shortest_path_matrix, get_random_walk_matrix, get_resistance

def find_connected_components(net):
    """
    Find all connected components in the network, filtering out small components.

    :param net: NetworkX graph.
    :return: List of components with more than 3 nodes.
    """
    return [list(comp) for comp in nx.connected_components(net) if len(comp) > 3]


def fisher_combine_pvalues(pvals: List[float]) -> float:
    p = np.clip(np.asarray(pvals, float), 1e-300, 1.0)  # avoid log(0)
    X = -2.0 * np.sum(np.log(p))
    return float(chi2.sf(X, 2 * p.size))


def qap_corr(A: np.ndarray, B: np.ndarray, perms: int, seed: int = None) -> Tuple[float, float]:
    A = np.asarray(A, float); B = np.asarray(B, float)
    n = A.shape[0]; assert A.shape == B.shape == (n, n)
    I = np.eye(n, dtype=bool)
    ii, jj = np.where(~I)  # off-diagonal indices

    y = B[ii, jj] - B[ii, jj].mean()
    y_norm = np.linalg.norm(y)
    if y_norm == 0:
        return np.nan, 1.0

    def _corr(x):
        x = x - x.mean()
        xn = np.linalg.norm(x)
        return 0.0 if xn == 0 else float(np.dot(x, y) / (xn * y_norm))

    r_obs = _corr(A[ii, jj])

    rng = np.random.default_rng(seed)
    ge = 0
    for _ in range(perms):
        p = rng.permutation(n)
        r_b = _corr(A[p[ii], p[jj]])
        # two-sided tail
        ge += (r_b >= abs(r_obs)) if r_obs >= 0 else (r_b <= -abs(r_obs))
    p_val = (1 + ge) / (1 + perms)
    return r_obs, p_val


def calculate_qap(net, fst_matrix: np.ndarray, dist_type: str, perms: int):
    if dist_type == 'euclidean':
        D = get_euclidean_matrix(net)
    elif dist_type == 'path':
        D = get_shortest_path_matrix(net)
    elif dist_type == 'random':
        D = get_mfpt(net)
    elif dist_type == 'resistance':
        D = get_resistance(net)
    else:
        raise ValueError(f"Unknown distance type: {dist_type}")

    if nx.is_connected(net):
        return qap_corr(D, fst_matrix, perms)

    r_vals, p_vals, w = [], [], []
    for comp_nodes in find_connected_components(net):
        r, p = qap_corr(D[np.ix_(comp_nodes, comp_nodes)],
                        fst_matrix[np.ix_(comp_nodes, comp_nodes)],
                        perms)
        r_vals.append(r)
        p_vals.append(p)
        w.append(len(comp_nodes))

    if not r_vals:
        return None

    r = float(np.average(r_vals, weights=w))
    p = fisher_combine_pvalues(p_vals)
    return r, p


def qap_step(step, net, fst, dist_type, perms, replica):
    """
    Process a single step for Mantel test.
    :param step: Step index
    :param net: Network graph for the step
    :param fst: FST matrix for the step
    :param dist_type: Type of distance metric ('euclidean', 'path', 'random', 'resistance')
    :param perms: Number of permutations for Mantel test
    :param replica: Replica index
    :return: Mantel correlation result
    """
    result = calculate_qap(net, fst, dist_type, perms)
    if result is None:
        return None
    r, p = result
    return {'step': step, 'r_val': r, 'p_val': p, 'replica': replica}


def calculate_qap_process(data, perms, dist_type, replica, num_workers=None):
    """
    Calculate Mantel correlation and p-value for each step along fragmentation.
    :param data: Fragmentation data
    :param perms: Number of permutations for Mantel test
    :param dist_type: Distance type ('euclidean', 'path', 'random', 'resistance')
    :param replica: Replica index
    :param num_workers: Number of threads for parallelization
    :return: DataFrame with Mantel correlation results
    """
    results = []
    networks = data.networks[replica]
    fst_matrices = data.fst_matrices[replica]

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(qap_step, step, net, fst, dist_type, perms, replica)
                   for step, (net, fst) in enumerate(zip(networks, fst_matrices))]

        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result is not None:
                results.append(result)

    return pd.DataFrame(results)


def calculate_qap_replicas(data, perms, dist_type):
    """
    Calculate Mantel correlation and p-value across fragmentation for all replicas.
    :param data: Data for all fragmentation types
    :param perms: Number of permutations for Mantel test
    :param dist_type: Distance metric type ('euclidean', 'path', 'random', 'resistance')
    :return: DataFrame with Mantel correlation results for all replicas
    """
    results = []
    networks = data.networks

    for replica in range(len(networks)):
        print(f"Processing replica {replica}")
        cor_data = calculate_qap_process(data, perms, dist_type, replica)
        results.append(cor_data)

    return pd.concat(results)


def calculate_qap_all(data, perms, dist_type: str):
    """
    Calculate Mantel correlation and p-value across fragmentation for all fragmentation types.
    :param data: Data for all fragmentation types
    :param perms: Number of permutations for Mantel test
    :param dist_type: Distance metric type ('euclidean', 'path', 'random', 'resistance')
    :return: DataFrame with Mantel correlation results for all fragmentation types
    """
    results = []
    for frag_type in data.keys():
        print(f"Processing fragmentation type: {frag_type}")
        cor_data = calculate_qap_replicas(data[frag_type], perms, dist_type)
        cor_data['fragmentation_type'] = frag_type
        results.append(cor_data)

    cor_data = pd.concat(results)

    cor_data.to_csv(f'fst_{dist_type}.csv', index=False)
    return cor_data
