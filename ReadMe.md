#  Network-based genetic monitoring of landscape fragmentation

The repository includes a Python-based toolkit developed for analyzing the genetic effects of changes in connectivity in population networks.
For more information read the full article of Peled et al. at DOI

## Sections

### genetic_metrics 
Pipeline for calculating and analyzing population genetic metrics (heterozygosity, F_ST) from network fragmentation simulations.

- **`run_pipeline.py`** – Entry point script that generates RGG networks, applies all 8 fragmentation types (rand, cor, intr, reg, div, dist, opt, wrst) in parallel, and serializes results to pickle files.

- **`funcs_run_pipeline.py`** – Core pipeline logic: maps fragmentation keys to functions, runs single-replicate fragmentation (graph → migration sequence → coalescence/F_ST matrices → het/fst distributions), and orchestrates parallel execution over replicates.

- **`mean_genetics.py`** – Post-processing utilities: normalizes and combines replicate-level heterozygosity and F_ST data across fragmentation types, and plots mean ± SD over fragmentation steps.

- **`Transformation.py`** – Wrapper around C library (`libmigration.so`) for transforming migration matrices into coalescence time and F_ST matrices using ctypes; defines `Migration` and `Coalescence` classes.

- **`libmigration.c`** – C implementation using GSL for solving linear systems: converts migration matrices to coefficient matrices and computes coalescence times.



### networks

Utilities for generating random-geometric graphs (RGGs) and computing network topology metrics across fragmentation steps.

- **`networks_generator.py`** – Creates random-geometric graphs with target edge counts and converts undirected networks to directed/asymmetric weighted graphs using mass-balanced sampling (conservative migration constraint: Σⱼ mᵢⱼ = Σⱼ mⱼᵢ).

- **`network_structures.py`** – Computes topological metrics (giant component fraction, isolated nodes, secondary components, "waste" components) per graph and aggregates over fragmentation steps and replicates; includes stacked-area plotting.


### fragmentation

Implements 8 different edge-removal strategies for simulating network fragmentation scenarios.

- **`processes.py`** – Defines all fragmentation algorithms, each returning a sequence of progressively fragmented networks:
  - **`remove_edge_random`** – Randomly delete bidirectional edges until none remain.
  - **`remove_edge_intrusive`** – Pick a random connected node; remove all its incident edges; repeat.
  - **`remove_edge_correlated`** – Always remove an edge adjacent to the previously deleted one; fall back to random when isolated.
  - **`remove_edge_distance`** – Remove edges from longest to shortest (requires `pos` node attribute).
  - **`remove_edge_regressive`** – Progressively strip edges starting with western-most nodes (requires `pos`).
  - **`remove_edge_divisive`** – Draw random border-to-border dividing lines; remove intersected edges west-to-east; repeat.
  - **`remove_edge_optimal`** – Re-compute edge betweenness centrality after each removal; drop the lowest-betweenness edge.
  - **`remove_edge_worst`** – Drop the edge with maximum betweenness centrality at each step.


### early_warning

Functions for analyzing heterozygosity trajectories at late fragmentation stages to identify early warning signals of population collapse.

- **`early_warning.py`** – Tools for detecting critical transitions in fragmented populations:

### correlations

Functions for testing and visualizing correlations between genetic distance (F_ST), network topology (centrality), and spatial/graph-theoretic distance metrics.

- **`run_cor.py`** – Entry point script: loads fragmentation data, computes QAP tests for F_ST vs. distance matrices (resistance/path/etc.), and generates plots.

- **`cor_test.py`** – Implements QAP (Quadratic Assignment Procedure) correlation tests between F_ST and network distance matrices (Euclidean, shortest path, random walk, resistance). Handles disconnected graphs by computing per-component correlations weighted by size and combining p-values via Fisher's method.

- **`cor_plot.py`** – Plotting utilities for Mantel correlation results: line plots of correlation vs. fragmentation step (filtered by p < 0.05), and scatter plots of F_ST vs. distance for individual steps.

- **`mantel_test.py`** – Mantel test implementation (wrapper around `mantel` package) for correlating F_ST and distance matrices, with support for disconnected networks via weighted component-wise aggregation. Includes parallel processing for multi-step/multi-replica analyses.

- **`distance_matrices.py`** – Computes pairwise distance matrices from networks:

- **`centrality.py`** – Computes node-level centrality measures (degree, betweenness) for single networks, across replicates, and across fragmentation types; exports to CSV.

- **`centrality_corr.py`** – Merges centrality data with heterozygosity data, computes Pearson correlations between centrality measures and heterozygosity per (frag_type, replica, step), and exports results.

- **`nodes_matrices.py`** – Utilities for node-level correlation analysis.

### other

Miscellaneous analysis and plotting utilities for exploring heterozygosity distributions, variance, giant component relationships, and robustness tests.

- **`distributions.py`** – Tools for visualizing heterozygosity/F_ST distributions across fragmentation steps:

- **`giant_comp.py`** – Analyzes relationship between giant component size and heterozygosity:

- **`pop_ind.py`** – Individual node trajectory analysis:

- **`robustness_tests.py`** – Robustness validation for conservative migration balancing (`project_to_conservative`): generates random fragmentation replicates, tracks migration matrix entries before/after balancing at each step, and exports detailed tidy DataFrame of changes.

- **`variance.py`** – Heterozygosity variance analysis:


## Contact
This repository is provided for academic and research purposes.
For questions, suggestions, or contributions, please contact ohad.peled@mail.huji.ac.il


