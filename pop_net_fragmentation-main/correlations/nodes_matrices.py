from typing import List, Dict, Literal
import networkx as nx
import numpy as np
import pandas as pd
from infomap import Infomap
from joypy import joyplot
from matplotlib import pyplot as plt
import seaborn as sns
from mantel import test
from scipy.stats import pearsonr
from scipy.stats import norm
from funcs import load_data, FragmentationResult, assign_node_numbers, percent_step







############ correlation of sigle steps
# def add_annotation(ax, r: float, p: float) -> None:
#     """
#     Add annotation with correlation coefficient and p-value to the plot.
#
#     :param ax: Matplotlib Axes object.
#     :param r: Pearson correlation coefficient.
#     :param p: P-value of the correlation.
#     """
#     annotation_text = f'r = {r:.2f}\np = {p:.3f}' if p >= 0.001 else f'r = {r:.2f}\np < 0.001'
#     ax.annotate(
#         annotation_text,
#         xy=(0.5, 0.1),  # Position (x, y) as relative plot coordinates
#         xycoords='axes fraction',  # Use axes fraction for relative positioning
#         fontsize=16,
#         style='italic',
#         fontname='serif'
#     )
# def plot_correlation(
#     df: pd.DataFrame,
#     measure: Literal['degree', 'betweenness'],
#     output_path: str
# ) -> None:
#     """
#     Plot the correlation between a centrality measure and heterozygosity for
#      fragmentation type-step-replica.
#
#     :param df: DataFrame containing the centrality and heterozygosity data.
#     :param measure: The centrality measure to correlate ('degree_centrality' or 'betweenness_centrality').
#     :param output_path: Path to save the plot.
#     """
#     # 1. Compute the correlation coefficient (Pearson)
#     r, p = pearsonr(df[measure], df['het'])
#
#     plt.figure(figsize=(6, 4))
#     sns.regplot(data=df, x=measure, y='het', fit_reg=True)
#
#     # 3. Annotate the plot with r and p-value
#     add_annotation(plt.gca(), r, p)
#
#     # 3. Add the correlation coefficient to the plot
#     plt.xlabel('Degree', fontsize=18)
#     plt.ylabel('Heterozygosity', fontsize=18)
#     plt.tick_params(axis='both', which='major', labelsize=14)
#     plt.ylim(-0.05, 1.2)
#
#     plt.savefig(output_path, format='svg')
#     plt.show()
# def preprocess_centrality_data(df: pd.DataFrame, replica: int, step: int, frag_type: str) -> pd.DataFrame:
#     """
#     Preprocess the centrality DataFrame by filtering for a specific replica, step, and fragmentation type.
#
#     :param df: DataFrame containing centrality data.
#     :param replica: Replica index to filter.
#     :param step: Step index to filter.
#     :param frag_type: Fragmentation type to filter.
#     :return: Filtered DataFrame.
#     """
#     filtered_df = df[(df['replica'] == replica) & (df['step'] == step)]
#     filtered_df = filtered_df[filtered_df['frag_type'] == frag_type]
#     return filtered_df
# def plot_correlation_steps(
#     df: pd.DataFrame,
#     frag_type: str,
#     replica: int,
#     steps: List[int],
#     measure: Literal['degree', 'betweenness'],
#     output_path: str
# ) -> None:
#     """
#     Produce a row of three scatter+regression plots of centrality vs. het,
#     for a single frag_type and replica, at the specified steps.
#
#     :param df: DataFrame with columns ['frag_type','replica','step','node_number',
#                'degree','betweenness','het'].
#     :param frag_type: Fragmentation type to filter on.
#     :param replica: Replica index to filter on.
#     :param steps: step indices to plot.
#     :param measure: Which centrality to plot ('degree' or 'betweenness').
#     :param output_path: Where to save the combined figure.
#     """
#     # set up 1×3 axes
#     fig, axes = plt.subplots(1, 3, figsize=(10, 2), sharey=True)
#     for ax, step in zip(axes, steps):
#         # filter for frag_type, replica, and step
#         sub = preprocess_centrality_data(df, replica, step, frag_type)
#
#         sns.regplot(
#             data=sub,
#             x=measure,
#             y='het',
#             ax=ax,
#             scatter_kws={'alpha':0.7},
#         )
#
#         # compute and annotate r & p
#         r, p = pearsonr(sub[measure], sub['het'])
#         add_annotation(ax, r, p)
#
#         # styling
#         ax.set_xlabel(measure.capitalize(), fontsize=14)
#         if ax is axes[0]:
#             ax.set_ylabel("Heterozygosity", fontsize=14)
#         else:
#             ax.set_ylabel("")
#
#         ax.tick_params(labelsize=12)
#         ax.set_ylim(-0.1, 1.4)
#
#     plt.savefig(output_path, format='svg')
#     plt.show()
#
##### plot centrality vs. heterozygosity
# centrality_df = pd.read_csv('./csv_new/centrality_het.csv')
# steps = [0, 75, 150]
# plot_correlation_steps(
#     df=centrality_df,
#     frag_type='dist',
#     replica=10,
#     steps=steps,
#     measure='betweenness',
#     output_path='./figs/het_bet_steps.svg'
# )
