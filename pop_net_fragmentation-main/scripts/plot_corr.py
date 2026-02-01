import pandas as pd

from calculations.node_matrices import compute_het_central_correlation
from data_manipulation.manp_node_matrices import merge_centrality_het, filter_correlations
from funcs import load_data
from viz_funcs.viz_node_matrices import plot_correlation

###scripts
##### compute centrality for all fragmentation types
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
data = load_data(fragmentation_types)


def compute_centralities_types(data, fragmentation_types):
	pass


centrality_df = compute_centralities_types(data, fragmentation_types)

#### merge centrality with heterozygosity data
frag_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
data = load_data(frag_types)
centrality_df = pd.read_csv('./csv_new/centrality.csv')
merged_df = merge_centrality_het(centrality_df, data, frag_types)

#### compute correlation between centrality and heterozygosity
centrality_df = pd.read_csv('./csv_new/centrality_het.csv')
corr_df = compute_het_central_correlation(
    df=centrality_df,
    centrality='betweenness',
)

#### plot correlation between centrality and heterozygosity
corr_df = pd.read_csv('./csv_new/het_bet_correlation.csv')
filtered_corr_df = filter_correlations(corr_df, min_replicates=5)
print(filtered_corr_df)
plot_correlation(
    corr_df=filtered_corr_df,
    output_path='./figs/corr_bet.svg'
)

