from mantel import test
from distance_matrices import get_euclidean_matrix, get_shortest_path_matrix, get_random_walk_matrix
import numpy as np
import networkx as nx
from typing import Optional
import concurrent.futures
import pandas as pd

def find_connected_components(net):
    """
    Find all connected components in the network, filtering out small components.

    :param net: NetworkX graph.
    :return: List of components with more than 3 nodes.
    """
    return [list(comp) for comp in nx.connected_components(net) if len(comp) > 3]




def calculate_mantel(net, fst_matrix, dist_type: str, perms: int) -> Optional[tuple]:
    """
    Perform a Mantel test for two distance matrices.
    :param net: NetworkX graph
    :param fst_matrix: Genetic distance matrix
    :param dist_type: Type of distance ('euclidean', 'path', or 'random')
    :param perms: Number of permutations for the Mantel test
    :return: Correlation coefficient (r) and p-value, or None if no valid correlation
    """
    # Select the appropriate distance matrix based on dist_type
    if dist_type == 'euclidean':
        distance_matrix = get_euclidean_matrix(net)
    elif dist_type == 'path':
        distance_matrix = get_shortest_path_matrix(net)
    elif dist_type == 'random':
        distance_matrix = get_random_walk_matrix(net)
    else:
        raise ValueError(f"Unknown distance type: {dist_type}")

    if nx.is_connected(net):
        r, p, _ = test(X=distance_matrix, Y=fst_matrix, perms=perms, method='pearson', ignore_nans=True)
        return r, p

    # Handle disconnected networks by computing Mantel for each component
    r_values, p_values, weights = [], [], []
    components = find_connected_components(net)

    for comp in components:
        comp_dist_matrix = distance_matrix[np.ix_(comp, comp)]
        comp_fst_matrix = fst_matrix[np.ix_(comp, comp)]
        r, p, _ = test(X=comp_dist_matrix, Y=comp_fst_matrix, perms=perms, method='pearson', ignore_nans=True)

        r_values.append(r)
        p_values.append(p)
        weights.append(len(comp))  # Weight by component size

    if not r_values or not p_values:
        return None

    # Calculate weighted averages of r and p across components
    weighted_r = np.average(r_values, weights=weights)
    weighted_p = np.average(p_values, weights=weights)

    return weighted_r, weighted_p


def mantel_step(step, net, fst, dist_type, perms, replica):
    """
    Process a single step for Mantel test.
    :param step: Step index
    :param net: Network graph for the step
    :param fst: FST matrix for the step
    :param dist_type: Type of distance metric ('euclidean', 'path', 'random')
    :param perms: Number of permutations for Mantel test
    :param replica: Replica index
    :return: Mantel correlation result
    """
    result = calculate_mantel(net, fst, dist_type, perms)
    if result is None:
        return None
    r, p = result
    return {'step': step, 'r_val': r, 'p_val': p, 'replica': replica}


def calculate_mantel_process(data, perms, dist_type, replica, num_workers=None):
    """
    Calculate Mantel correlation and p-value for each step along fragmentation.
    :param data: Fragmentation data
    :param perms: Number of permutations for Mantel test
    :param dist_type: Distance type ('euclidean', 'path', 'random')
    :param replica: Replica index
    :param num_workers: Number of threads for parallelization
    :return: DataFrame with Mantel correlation results
    """
    results = []
    networks = data.networks[replica]
    fst_matrices = data.fst_matrices[replica]

    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(mantel_step, step, net, fst, dist_type, perms, replica)
                   for step, (net, fst) in enumerate(zip(networks, fst_matrices))]

        for future in concurrent.futures.as_completed(futures):
            result = future.result()
            if result is not None:
                results.append(result)

    return pd.DataFrame(results)


def calculate_mantel_replicas(data, perms, dist_type):
    """
    Calculate Mantel correlation and p-value across fragmentation for all replicas.
    :param data: Data for all fragmentation types
    :param perms: Number of permutations for Mantel test
    :param dist_type: Distance metric type ('euclidean', 'path', 'random')
    :return: DataFrame with Mantel correlation results for all replicas
    """
    results = []
    networks = data.networks

    for replica in range(len(networks)):
        print(f"Processing replica {replica}")
        cor_data = calculate_mantel_process(data, perms, dist_type, replica)
        results.append(cor_data)

    return pd.concat(results)


def calculate_mantel_all(data, perms, dist_type: str):
    """
    Calculate Mantel correlation and p-value across fragmentation for all fragmentation types.
    :param data: Data for all fragmentation types
    :param perms: Number of permutations for Mantel test
    :param dist_type: Distance metric type ('euclidean', 'path', 'random')
    :return: DataFrame with Mantel correlation results for all fragmentation types
    """
    results = []
    for frag_type in data.keys():
        print(f"Processing fragmentation type: {frag_type}")
        cor_data = calculate_mantel_replicas(data[frag_type], perms, dist_type)
        cor_data['fragmentation_type'] = frag_type
        results.append(cor_data)

    cor_data = pd.concat(results)

    cor_data.to_csv(f'./csv_new/fst_{dist_type}_corrrlation.csv', index=False)
    return cor_data
