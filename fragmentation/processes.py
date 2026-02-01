from __future__ import annotations

import random
from math import sqrt
from typing import Iterable, List, Set, Tuple

import networkx as nx

__all__ = [
    "remove_edge_random",
    "remove_edge_intrusive",
    "remove_edge_correlated",
    "remove_edge_distance",
    "remove_edge_regressive",
    "remove_edge_divisive",
    "remove_edge_optimal",
    "remove_edge_worst",
]

# --------------------------------------------------------------------------- #
# Utility helpers
# --------------------------------------------------------------------------- #
RNG = random.Random()  # module‑level PRNG (override via seed args when needed)


def _remove_edge_pair(net: nx.Graph, u: int, v: int) -> None:
    """Delete (u, v) and (v, u) if present."""
    if net.has_edge(u, v):
        net.remove_edge(u, v)
    if net.is_directed() and net.has_edge(v, u):
        net.remove_edge(v, u)


def _connected_nodes(net: nx.Graph) -> Set[int]:
    """Nodes that belong to components of size > 1."""
    comps: Iterable[Set[int]] = (
        nx.weakly_connected_components(net)
        if net.is_directed()
        else nx.connected_components(net)
    )
    return {n for comp in comps if len(comp) > 1 for n in comp}


def _connected_edges(net: nx.Graph, nodes: Set[int]) -> List[Tuple[int, int]]:
    """Edges touching at least one node in *nodes*."""
    return [(u, v) for u, v in net.edges() if u in nodes or v in nodes]


# ---- geometry helper for divisive removal ------------------------------------

def _segments_intersect(p: Tuple[float, float], q: Tuple[float, float], r: Tuple[float, float], s: Tuple[float, float]) -> bool:
    """Return True if segments pq and rs intersect (collinear counts)."""

    def orient(a, b, c):
        return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])

    def on_segment(a, b, c):
        return min(a[0], b[0]) <= c[0] <= max(a[0], b[0]) and min(a[1], b[1]) <= c[1] <= max(a[1], b[1])

    o1, o2, o3, o4 = orient(p, q, r), orient(p, q, s), orient(r, s, p), orient(r, s, q)
    if o1 * o2 < 0 and o3 * o4 < 0:
        return True
    return (
        (o1 == 0 and on_segment(p, q, r))
        or (o2 == 0 and on_segment(p, q, s))
        or (o3 == 0 and on_segment(r, s, p))
        or (o4 == 0 and on_segment(r, s, q))
    )


def _intersection_x(p: Tuple[float, float], q: Tuple[float, float], r: Tuple[float, float], s: Tuple[float, float]) -> float:
    """x‑coordinate where lines pq and rs intersect (assumes they do)."""
    px, py, qx, qy = *p, *q
    rx, ry, sx, sy = *r, *s
    denom = (qx - px) * (sy - ry) - (qy - py) * (sx - rx)
    if abs(denom) < 1e-12:  # collinear – choose min x of the segment as proxy
        return min(rx, sx)
    t = ((rx - px) * (sy - ry) - (ry - py) * (sx - rx)) / denom
    return px + t * (qx - px)


# --------------------------------------------------------------------------- #
# Fragmentation scenarios
# --------------------------------------------------------------------------- #

def remove_edge_random(G: nx.Graph, *, seed: int | None = None) -> List[nx.Graph]:
    """Randomly delete bidirectional edges until none remain."""
    rnd = random.Random(seed) if seed is not None else RNG
    net = G.copy()
    net_list = [net.copy()]

    while net.number_of_edges():
        u, v = rnd.choice(list(net.edges()))
        _remove_edge_pair(net, u, v)
        net_list.append(net.copy())

    return net_list


def remove_edge_intrusive(G: nx.Graph, *, seed: int | None = None) -> List[nx.Graph]:
    """Pick a random connected node; remove all its incident edges; repeat."""
    rnd = random.Random(seed) if seed is not None else RNG
    net = G.copy()
    net_list = [net.copy()]

    while net.number_of_edges():
        active = [n for n, deg in net.degree() if deg]
        if not active:
            break
        node = rnd.choice(active)
        for _, nbr in list(net.edges(node)):
            _remove_edge_pair(net, node, nbr)
            net_list.append(net.copy())

    return net_list


def remove_edge_correlated(G: nx.Graph, *, seed: int | None = None) -> List[nx.Graph]:
    """Always remove an edge adjacent to the previously deleted one; fall back to a random edge when isolated."""
    rnd = random.Random(seed) if seed is not None else RNG
    net = G.copy()
    net_list = [net.copy()]

    if not net.edges():
        return net_list

    u, v = rnd.choice(list(net.edges()))
    while net.number_of_edges():
        candidates = list(net.edges(u)) + list(net.edges(v))
        if not candidates:
            nodes = _connected_nodes(net)
            candidates = _connected_edges(net, nodes)
            if not candidates:
                break
        u, v = rnd.choice(candidates)
        _remove_edge_pair(net, u, v)
        net_list.append(net.copy())

    return net_list


def remove_edge_distance(G: nx.Graph) -> List[nx.Graph]:
    """Remove edges from longest to shortest (requires 'pos' node attribute)."""
    net = G.copy()
    pos = nx.get_node_attributes(net, "pos")
    net_list = [net.copy()]

    def length(edge: Tuple[int, int]) -> float:
        (x1, y1), (x2, y2) = pos[edge[0]], pos[edge[1]]
        return sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

    for u, v in sorted(net.edges(), key=length, reverse=True):
        _remove_edge_pair(net, u, v)
        yield_net = net.copy()
        # building list incrementally avoids extra pass over edges
        net_list.append(yield_net)

    return net_list


def remove_edge_regressive(G: nx.Graph) -> List[nx.Graph]:
    """Progressively strip edges starting with western‑most nodes (needs 'pos')."""
    net = G.copy()
    pos = nx.get_node_attributes(net, "pos")
    x_order = sorted(pos, key=lambda n: pos[n][0])
    net_list = [net.copy()]

    while net.number_of_edges():
        for n in x_order:
            edges = list(net.edges(n))
            if edges:
                edges.sort(key=lambda e: x_order.index(e[1] if e[0] == n else e[0]))
                u, v = edges[0]
                _remove_edge_pair(net, u, v)
                net_list.append(net.copy())
                break
        else:
            break
    return net_list

def _border_divider(rnd: random.Random) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    """Return two points such that the segment spans the unit square borders."""
    while True:
        # pick one point on a vertical border, one point on a horizontal border
        x1 = rnd.choice([0.0, 1.0])
        y1 = rnd.random()
        x2 = rnd.random()
        y2 = rnd.choice([0.0, 1.0])
        if abs(x1 - x2) > 1e-9 and abs(y1 - y2) > 1e-9:  # avoid pure vertical / horizontal
            return (x1, y1), (x2, y2)


def remove_edge_divisive(G: nx.Graph, *, seed: int | None = None) -> List[nx.Graph]:
    rnd = random.Random(seed) if seed is not None else RNG
    net = G.copy()
    pos = nx.get_node_attributes(net, "pos")
    net_list = [net.copy()]

    while net.number_of_edges():
        p1, p2 = _border_divider(rnd)
        affected = [
            e for e in net.edges()
            if _segments_intersect(p1, p2, pos[e[0]], pos[e[1]])
        ]
        if not affected:
            continue  # try a new divider
        # Sort by actual intersection x‑coordinate (west → east)
        affected.sort(key=lambda e: _intersection_x(p1, p2, pos[e[0]], pos[e[1]]))
        for u, v in affected:
            _remove_edge_pair(net, u, v)
            net_list.append(net.copy())
    return net_list


def remove_edge_optimal(G: nx.Graph) -> List[nx.Graph]:
    """Re‑evaluate edge betweenness after every removal; drop the lowest‑betweenness edge each step."""
    net = G.copy()
    net_list = [net.copy()]

    while net.number_of_edges():
        centrality = nx.edge_betweenness_centrality(net)
        u, v = min(centrality, key=centrality.get)
        _remove_edge_pair(net, u, v)
        net_list.append(net.copy())

    return net_list


def remove_edge_worst(G: nx.Graph) -> List[nx.Graph]:
    """At each step drop the edge with *maximum* betweenness centrality."""
    net = G.copy()
    net_list = [net.copy()]

    while net.number_of_edges():
        (u, v), _ = max(nx.edge_betweenness_centrality(net).items(), key=lambda x: x[1])
        _remove_edge_pair(net, u, v)
        net_list.append(net.copy())

    return net_list
