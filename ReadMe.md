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


**Note:** The directory also imports from `stats_utils.py`, `processes.py`, `networks_generator.py`, and `funcs.py` located at the repository root. [View all files in this directory on GitHub](https://github.com/Greenbaum-Lab/pop_net_fragmentation/tree/main/genetic_metrics).
- **`Transformation.py`**  
  Implements mathematical transformations between migration matrices, coalescence times, and Fst statistics, including routines for component detection and conservative migration matrix generation, using both Python and optional C acceleration.

- **`processes.py`**  
  Contains core simulation logic for different fragmentation processes, including random, correlated, distance-based, and optimal edge removal routines.

- **`funcs_initial_data.py`**  
  Generates and normalizes initial networks, runs fragmentation replicates, and summarizes genetic statistics for various network models and fragmentation types.

- **`run_pipeline.py`**  
  Script to run the full fragmentation analysis pipeline, including network creation, simulation for all fragmentation types, and batch pickle output.

- **`funcs.py`**  
  Defines data structures and utility functions for loading, processing, and summarizing fragmentation simulation results, including node/step annotation and component analysis.

- **`mean_genetics.py`**  
  Calculates and visualizes mean heterozygosity and Fst across fragmentation types, with normalization and grouped analysis.

- **`distributions.py`**  
  Extracts, filters, and visualizes distributions of genetic diversity measures (e.g., heterozygosity, Fst) at fixed fragmentation intervals.

- **`centrality.py`**  
  Computes node centralities (degree, betweenness) for single or multiple networks, aiding network fragmentation analysis.

- **`centrality_corr.py`**  
  Analyzes correlations between node centrality measures and heterozygosity, with tools for merging, filtering, and plotting these relationships across fragmentation experiments.

- **`distance_matrices.py`**  
  Calculates shortest path, Euclidean, and random walk distance matrices for network nodes using NetworkX and parallel processing.

- **`giant_comp.py`**  
  Analyzes the relationship between the size of the largest network component and genetic diversity, including binning and plotting tools.

- **`mantel_plot.py`**  
  Plots Mantel test results and correlations between genetic and distance matrices across fragmentation types and steps.

- **`mantel_test.py`**  
  Performs Mantel tests to assess the correlation between genetic and distance matrices, supporting parallel and type-specific analyses.

- **`network_structures.py`**  
  Computes, aggregates, and visualizes network structure metrics (giant component, isolated nodes, secondary components, waste) for single and multiple replicates.

- **`nodes_matrices.py`**  
  Provides tools for analyzing node-level statistics and correlations, including centrality and heterozygosity relationships.

- **`pop_ind.py`**  
  Utilities for selecting and plotting heterozygosity trajectories of individual nodes in population fragmentation simulations.

- **`variance`**  
  Computes and plots the variance of heterozygosity across fragmentation types and steps, returning per-replica and aggregated variance statistics.

- **`early_warning.py`**  
  Computes early warning indicators for genetic collapse or transitions (e.g., standard deviation, skewness, kurtosis, return rate) and provides plotting utilities.

- **`libmigration.c` / `libmigration.so`**  
  Supplies C routines and a compiled library for efficient migration, coalescence, and Fst matrix computations, used as an optional backend for high-performance operations.



## Contact
This repository is provided for academic and research purposes.
For questions, suggestions, or contributions, please contact ohad.peled@mail.huji.ac.il

