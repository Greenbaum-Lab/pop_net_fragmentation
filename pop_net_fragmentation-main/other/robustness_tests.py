# """
# Generate 10 replicas of the RANDOM fragmentation process and track how
# project_to_conservative (PC) adjusts every migration entry.
#
# Output file
# -----------
# rand_matrix_pre_post_REP10_TEST.pickle  containing a tuple:
#     • all_pre   – list[rep] of list[step] of np.ndarray (M_pre)
#     • all_post  – list[rep] of list[step] of np.ndarray (M_post)
#     • df_pairs  – pandas DataFrame with columns
#           replicate, step, source, target,
#           w_before, w_after, delta_pc
# """
#
# import sys, pickle, random, pathlib
# from copy import deepcopy
#
# import numpy  as np
# import pandas as pd
# import networkx as nx
#
# # ───────────── project modules ──────────────────────────────────────────
#
#
# from networks_generator import make_rgg, project_to_conservative
# from processes          import _remove_edge_pair        # bidirectional deletion
#
# # ───────────── helper ───────────────────────────────────────────────────
# def to_matrix(g: nx.DiGraph) -> np.ndarray:
#     """Return full n×n numpy migration matrix (nodes sorted by index)."""
#     return nx.to_numpy_array(g, dtype=float, weight="weight",
#                              nodelist=sorted(g.nodes()))
#
# # ───────────── parameters ───────────────────────────────────────────────
# n_replicates   = 10
# n_nodes        = 50
# target_edges   = 250
# base_seed      = 42      # will offset by +rep to keep replicates independent
#
# # ───────────── main loop ────────────────────────────────────────────────
# all_pre, all_post, df_chunks = [], [], []
#
# for rep in range(n_replicates):
#     # separate NumPy and Python RNGs for reproducibility
#     rng_np = np.random.default_rng(base_seed + rep)
#     rng_py = random.Random(base_seed + rep)
#
#     # 1. initial asymmetric network (already directed)
#     G = make_rgg(1, n_nodes, target_edges,
#                  asymmetric=True, rng=rng_np)[0]
#     project_to_conservative(G)
#
#     net = deepcopy(G)
#     mats_pre, mats_post = [], []
#     step = 0
#
#     while net.number_of_edges():
#         # ---- edge removal -------------------------------------------------
#         u, v = rng_py.choice(list(net.edges()))
#         _remove_edge_pair(net, u, v)
#
#         M_pre  = to_matrix(net)          # after removal
#         mats_pre.append(M_pre.copy())
#
#         # ---- balance ------------------------------------------------------
#         project_to_conservative(net)
#         M_post = to_matrix(net)          # after PC
#         mats_post.append(M_post.copy())
#
#         # ---- flatten to tidy DataFrame chunk ------------------------------
#         n = M_pre.shape[0]
#         src, tgt = np.meshgrid(range(n), range(n), indexing="ij")
#         df_chunks.append(pd.DataFrame({
#             "replicate": rep,
#             "step":      step,
#             "source":    src.ravel(),
#             "target":    tgt.ravel(),
#             "w_before":  M_pre.ravel(),
#             "w_after":   M_post.ravel(),
#         }))
#         step += 1
#
#     all_pre.append(mats_pre)
#     all_post.append(mats_post)
#
# # full tidy table with PC deltas
# df_pairs = pd.concat(df_chunks, ignore_index=True)
# df_pairs["delta_pc"] = df_pairs["w_after"] - df_pairs["w_before"]
#
# # ───────────── save everything ──────────────────────────────────────────
# out_file = "rand_matrix_pre_post_REP10_TEST.pickle"
# with open(out_file, "wb") as f:
#     pickle.dump((all_pre, all_post, df_pairs), f)
#
# print(f"✔  Generated {n_replicates} replicas "
#       f"({df_pairs['replicate'].nunique()} unique) and saved to '{out_file}'.")
#
#
#
#
#
#
#
# # -------------------------------------------------------------
# # 0.  Load data
# # -------------------------------------------------------------
# import pickle, numpy as np, pandas as pd
# import matplotlib.pyplot as plt
#
# PKL = "rand_matrix_pre_post_REP10_TEST.pickle"   # output from previous script
#
# with open(PKL, "rb") as f:
#     all_pre, all_post, df_pairs = pickle.load(f)
#
# n_reps   = len(all_pre)
# n_nodes  = all_pre[0][0].shape[0]
# rng_ctrl = np.random.default_rng(999)
#
# # -------------------------------------------------------------
# # 1.  Build RANDOM-reassignment deltas for every (rep, step)
# # -------------------------------------------------------------
# rand_records = []
#
# for rep, mats_pre in enumerate(all_pre):
#     for step, M_pre in enumerate(mats_pre):
#         rand_w  = np.abs(rng_ctrl.normal(1.0, 0.3, size=M_pre.shape))
#         M_rand  = np.where(M_pre > 0, rand_w, 0.0)
#
#         mask = M_pre > 0                          # edges still alive
#         mean_abs = np.abs(M_rand - M_pre)[mask].mean()
#         rand_records.append({"replicate": rep,
#                              "step":      step,
#                              "mean_abs_rand": mean_abs})
#
# df_rand_mean = pd.DataFrame(rand_records)
#
# # -------------------------------------------------------------
# # 2.  Mean |Δ| from project_to_conservative (already in df_pairs)
# # -------------------------------------------------------------
# df_pairs["abs_delta_pc"] = df_pairs["delta_pc"].abs()
#
# df_pc_mean = (df_pairs
#               .groupby(["replicate", "step"])
#               ["abs_delta_pc"]
#               .mean()
#               .reset_index()
#               .rename(columns={"abs_delta_pc": "mean_abs_pc"}))
#
# # merge PC & random means
# df_means = (df_pc_mean
#             .merge(df_rand_mean, on=["replicate", "step"]))
#
# # -------------------------------------------------------------
# # 3.  Grand mean ± SD across replicas  (per step)
# # -------------------------------------------------------------
# agg = (df_means
#        .groupby("step")
#        .agg(mean_pc   = ("mean_abs_pc",   "mean"),
#             sd_pc     = ("mean_abs_pc",   "std"),
#             mean_rand = ("mean_abs_rand", "mean"),
#             sd_rand   = ("mean_abs_rand", "std"))
#        .reset_index())
#
# # -------------------------------------------------------------
# # 4-A.  Plot: mean ± SD  (PC vs. random)
# # -------------------------------------------------------------
# plt.figure(figsize=(8,4))
#
# # project_to_conservative
# plt.plot(agg["step"], agg["mean_pc"], label="Corrected")
# plt.fill_between(agg["step"],
#                  agg["mean_pc"]-agg["sd_pc"],
#                  agg["mean_pc"]+agg["sd_pc"],
#                  alpha=.25)
#
# # random reassignment
# plt.plot(agg["step"], agg["mean_rand"], label="Random")
# plt.fill_between(agg["step"],
#                  agg["mean_rand"]-agg["sd_rand"],
#                  agg["mean_rand"]+agg["sd_rand"],
#                  alpha=.25)
#
# plt.xlabel("Fragmentation step")
# plt.ylabel("Mean |Δ| per edge")
# plt.title("Mean magnitude of edge-weight change for 10 replicas")
# plt.legend()
# plt.tight_layout()
# plt.show()
#
# # -------------------------------------------------------------
# # 4-B.  Histograms: pooled ∣Δ∣ distributions
# # -------------------------------------------------------------
# # pool across all replicas + steps
# abs_rand_all = pd.concat([
#     pd.Series(np.abs(
#         np.where(M_pre > 0, np.abs(rng_ctrl.normal(1, 0.3, size=M_pre.shape)), 0.0) - M_pre
#     ).ravel())
#     for mats_pre in all_pre
#     for M_pre in mats_pre
# ], ignore_index=True)
#
# abs_pc_all   = df_pairs["abs_delta_pc"].loc[lambda x: x > 0]
# abs_rand_all = abs_rand_all.loc[lambda x: x > 0]   # from previous pooling
# import seaborn as sns
#
# plt.figure(figsize=(6,4))
# sns.kdeplot(abs_pc_all, fill=True, alpha=0.5, label="Corrected", bw_adjust=1.2)
# sns.kdeplot(abs_rand_all, fill=True, alpha=0.5, label="Random", bw_adjust=1.2)
# plt.ylim(0, 10)
# plt.xlabel("|Δ|")
# plt.ylabel("Density")
# plt.title("Edge-weight change distributions for 10 replicas")
# plt.legend()
# plt.tight_layout()
# plt.show()
#
# #print the range and mean of the distributions
# print(f"Corrected: range = [{abs_pc_all.min():.3f}, {abs_pc_all.max():.3f}], mean = {abs_pc_all.mean():.3f}")
# print(f"Random:    range = [{abs_rand_all.min():.3f}, {abs_rand_all.max():.3f}], mean = {abs_rand_all.mean():.3f}")






# -------------------------------------------------------------------
# 0.  Parameters and imports
# -------------------------------------------------------------------
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from networks_generator import make_rgg
from Transformation      import conservative_from_normal

n_reps        = 20          # <-- change for more / fewer robustness runs
n_nodes       = 50
target_edges  = 250
mu, sigma     = 1.0, 0.3
lower, upper  = 0.1, 4.0

# -------------------------------------------------------------------
# 1.  Collect diagnostics across replicas
# -------------------------------------------------------------------
pool_imb_before, pool_imb_after, pool_abs_delta = [], [], []
max_imb_before = []
pool_w_cons = []        # balanced (non-zero) migration weights

for rep in range(n_reps):
    seed       = 1000 + rep
    rng        = np.random.default_rng(seed)

    # ---- binary support --------------------------------------------------
    G = make_rgg(1, n_nodes, target_edges, asymmetric=False, rng=rng)[0]
    support = nx.to_numpy_array(G, dtype=int)

    # ---- raw random draw (exactly the routine's first step) --------------
    x_raw = np.clip(rng.normal(mu, sigma, size=support.sum()), lower, upper)
    M_raw = np.zeros_like(support, dtype=float)
    idx   = np.argwhere(support)
    M_raw[idx[:,0], idx[:,1]] = x_raw

    # ---- balanced matrix -------------------------------------------------
    M_cons = conservative_from_normal(support, mu, sigma, lower, upper, seed=seed)
    pool_w_cons.extend(M_cons[support.astype(bool)].ravel())

    # ---- store diagnostics ----------------------------------------------
    pool_imb_before.extend(np.abs(M_raw.sum(1) - M_raw.sum(0)))
    pool_imb_after .extend(np.abs(M_cons.sum(1) - M_cons.sum(0)))

    delta = M_cons - M_raw
    pool_abs_delta.extend(np.abs(delta[support.astype(bool)]))

    max_imb_before.append(max(pool_imb_before[-n_nodes:]))  # last replica slice

# -------------------------------------------------------------------
# 2.  Plots
# -------------------------------------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(15,4))


# (A) imbalance distribution
# ---------------------------------------------------------------
#  A.  imbalance density  (zeros removed, log-x)
# ---------------------------------------------------------------
thr = 1e-12                            # treat tiny numbers as zero
imb_before_nz = np.array(pool_imb_before)[np.array(pool_imb_before) > thr]

# log-spaced bins from min to max
lo, hi = imb_before_nz.min(), imb_before_nz.max()
bins = np.logspace(np.log10(lo), np.log10(hi), num=40)

ax[0].hist(imb_before_nz, bins=bins, alpha=.8)
# ax[0].set_xscale("log")
ax[0].set_xlabel("|row_sum – col_sum|")
ax[0].set_ylabel("Count of populations")
ax[0].set_title("Population-level imbalance BEFORE balancing")


# (B) per-edge correction magnitude
ax[1].hist(pool_abs_delta, bins=60, density=True, alpha=.75, color="tab:purple")
ax[1].set_xlabel("|Δ| per edge")
ax[1].set_ylabel("Density")
ax[1].set_title("Magnitude of corrections")

# (C) replicate-wise max imbalance (before vs after)
ax[2].bar(range(n_reps), max_imb_before, alpha=.7, label="before")
ax[2].bar(range(n_reps), [0]*n_reps,   alpha=.7, label="after")
ax[2].set_xlabel("replicate")
ax[2].set_ylabel("max |imbalance|")
ax[2].set_title("Worst imbalance per replicate")
ax[2].legend()
plt.tight_layout()
plt.show()

#plot the distribution of values after balancing for all edges and replicates
w_cons_nz = [w for w in pool_w_cons if w > 1e-12]   # drop exact zeros
plt.figure(figsize=(6,4))
plt.hist(w_cons_nz, bins=60, density=True, alpha=.75, color="tab:green")
plt.xlabel("Balanced migration weight")
plt.ylabel("Density")
plt.title("Distribution of non-zero weights after balancing")
plt.tight_layout()
plt.show()