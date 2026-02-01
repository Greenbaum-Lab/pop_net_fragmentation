import concurrent
from itertools import combinations
from random import random

import networkx as nx
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

from funcs import load_data
import random
import concurrent.futures
import random
import numpy as np
import networkx as nx
from concurrent.futures import ThreadPoolExecutor

##fst-distance relationship

def get_shortest_path_matrix(net):
    """
    Calculate the shortest path length between all pairs of nodes using Floyd-Warshall.
    Unconnected nodes are marked with inf.

    :param net: NetworkX graph.
    :return: Distance matrix of edges between nodes.
    """
    n = nx.number_of_nodes(net)
    distance_matrix = np.full((n, n), np.inf)
    np.fill_diagonal(distance_matrix, 0)  # Distance to self is 0

    # Use Floyd-Warshall algorithm for all pairs shortest path
    path_lengths = nx.floyd_warshall_numpy(net)

    # Assign the result to the distance matrix
    distance_matrix[:] = path_lengths
    return distance_matrix


def get_euclidean_matrix(net):
    """
    Calculate the Euclidean distance between all pairs of nodes in the network.
    Unconnected nodes are marked with inf.

    :param net: NetworkX graph.
    :return: Distance matrix of Euclidean distances between nodes.
    """
    n = nx.number_of_nodes(net)
    positions = nx.get_node_attributes(net, 'pos')

    # Extract node positions into a NumPy array
    pos_array = np.array([positions[node] for node in sorted(net.nodes())])

    # Use broadcasting to compute pairwise Euclidean distances
    diff = pos_array[:, np.newaxis, :] - pos_array[np.newaxis, :, :]
    distance_matrix = np.sqrt(np.sum(diff ** 2, axis=-1))

    # Handle unconnected nodes by setting distances to inf
    if not nx.is_connected(net):
        for u, v in combinations(range(n), 2):
            if not nx.has_path(net, u, v):
                distance_matrix[u, v] = np.inf
                distance_matrix[v, u] = np.inf

    return distance_matrix


def random_walk(net, start, end, max_steps=100):
    """
    Perform a random walk between two nodes until the target is reached or max_steps is exceeded.

    :param net: NetworkX graph.
    :param start: Start node.
    :param end: End node.
    :param max_steps: Maximum number of steps before stopping.
    :return: Number of steps in the walk.
    """
    current_node = start
    steps = 0
    while current_node != end and steps < max_steps:
        neighbors = list(net.neighbors(current_node))
        if neighbors:
            current_node = random.choice(neighbors)
        steps += 1
    return steps


def compute_random_walk_distance(net, u, v, num_walks=50):
    """
    Compute the random walk distance between nodes u and v.

    :param net: NetworkX graph.
    :param u: Start node.
    :param v: End node.
    :param num_walks: Number of random walks to simulate.
    :return: Mean number of steps over all walks, or np.inf if no path exists.
    """
    if not nx.has_path(net, u, v):
        return u, v, np.inf  # No valid path, return np.inf as the distance

    # Perform random walks if a path exists
    with ThreadPoolExecutor() as executor:
        steps = list(executor.map(lambda _: random_walk(net, u, v), range(num_walks)))

    return u, v, np.mean(steps)


def get_random_walk_matrix(net, num_workers=None):
    """
    Calculate the random walk distance between all pairs of nodes in the network.

    :param net: NetworkX graph.
    :param num_workers: Number of threads to use for parallel processing.
    :return: Distance matrix of random walk distances.
    """
    n = nx.number_of_nodes(net)
    distance_matrix = np.full((n, n), np.inf)
    np.fill_diagonal(distance_matrix, 0)

    pairs = [(u, v) for u in range(n) for v in range(u + 1, n)]

    with ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(compute_random_walk_distance, net, u, v) for u, v in pairs]

        for future in futures:
            u, v, res = future.result()
            distance_matrix[u, v] = res
            distance_matrix[v, u] = res

    return distance_matrix





