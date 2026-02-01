import os

from distributions import filter_intervals, plot_distribution
from funcs import load_data
from mean_genetics import plot_genetics
from pop_ind import select_random_nodes, extract_nodes, plot_het_nodes

directory = "/home/lab2/PycharmProjects/fragmentation"
os.chdir(directory)
print(f"Current working directory: {os.getcwd()}")

def plot_genetic_data():
    """
    Plot genetic data for different fragmentation types.
    choose fragmentation types and measure (het or fst)
    """
    fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
    data = load_data(fragmentation_types)
    plot_genetics(data, measure= 'fst')

def plot_distributions():
    """
    Plot distributions for a single fragmentation type.
    choose fragmentation types and measure (het or fst)
    """
    fragmentation_types = ['rand']
    data = load_data(fragmentation_types)
    df = filter_intervals(data['rand'], measure='het', interval_pct=25)
    plot_distribution(df, measure='het', frag_type='rand')

def plot_individual_nodes():
    """
    Plot individual nodes for a single fragmentation type.
    choose fragmentation type
    """
    fragmentation_types = ['rand']
    data = load_data(fragmentation_types)
    df = data.get('rand').het_dist
    selected_nodes = select_random_nodes(df, per_replica=10)
    df_selected = extract_nodes(df, selected_nodes)
    plot_het_nodes(df_selected, n_nodes=10)

def main():
    plot_genetic_data()
    plot_distributions()
    plot_individual_nodes()

if __name__ == "__main__":
    main()



