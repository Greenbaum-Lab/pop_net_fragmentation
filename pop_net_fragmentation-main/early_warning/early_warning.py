import networkx as nx
import pandas as pd
from matplotlib import pyplot as plt
import random

from funcs import assign_node_numbers


def get_focal_step(df, num_nodes=10):
    """
    Gets the maximum step in each replica where there are at least 10 'het' values greater than 0.02.

    :param df: DataFrame with 'replica', 'step', and 'het' columns.
    :return: Dictionary with replicas as keys and maximum steps as values.
    """
    # Group by 'replica' and 'step', filter 'het' values greater than 0.02, and count them
    het_counts = df[df['het'] > 0.02].groupby(['replica', 'step']).size()

    # Filter the groups where the count is at least 10
    het_counts = het_counts[het_counts >= num_nodes]

    # For each replica, get the maximum step that satisfies the condition
    max_steps = het_counts.reset_index(level='step').groupby('replica')['step'].max()

    return max_steps.to_dict()


def get_max_het_nodes(df, num_nodes=10):
    """
    Gets the 10 nodes in each replica that have the highest 'het' in the 10th last step of each replica.

    :param df: DataFrame with 'replica', 'step', 'het', and 'node_number' columns.
    :return: Dictionary with replicas as keys and lists of node numbers with the highest 'het' in the 10th last step as values.
    """
    # Get the last step of each replica with n surviving nodes
    top_nodes = {}
    focal_step = get_focal_step(df, num_nodes=num_nodes)

    for replica, step in focal_step.items():
        # Select only the rows that belong to the 10th last step of the current replica
        df_step = df[(df['replica'] == replica) & (df['step'] == step)]

        # Sort the DataFrame by 'het' in descending order and get the 'node_number' of the top 10 rows
        top_nodes_replica = df_step.sort_values(by='het', ascending=False)['node_number'].head(num_nodes).tolist()
        top_nodes[replica] = top_nodes_replica

    return top_nodes


def extract_steps_for_nodes(df, max_het_nodes):
    """
    Extracts all steps for each node in its corresponding replica from the DataFrame.

    :param df: DataFrame with 'replica', 'step', 'het', and 'node_number' columns.
    :param max_het_nodes: Dictionary with replicas as keys and node numbers with the highest 'het' in the last step as values.
    :return: DataFrame with the extracted rows.
    """

    extracted_rows = pd.DataFrame()
    for replica, all_nodes in max_het_nodes.items():

        for node_number in all_nodes:
            df_replica_node = df[(df['replica'] == replica) & (df['node_number'] == node_number)]
            extracted_rows = pd.concat([extracted_rows, df_replica_node], ignore_index=True)
    return extracted_rows


def export_het_csv(data, frag: str):
    """
    Export the heterozygosity data to a CSV file.
    :param data: full fragmentation data
    :param frag: fragmentation type
    :return: csv file
    """
    data = assign_node_numbers(data)
    surviving_nodes = get_max_het_nodes(data, num_nodes=10)
    final_df = extract_steps_for_nodes(data, surviving_nodes)
    final_df.to_csv(f'{frag}_het.csv')
    print('File saved successfully')
    return final_df


def get_largest_component(nets):
    """
    get the largest component of each network in each replica and step.
    for early warning analysis.
    :param nets: list of lists of networks in the format of nets[replica][step]
    :return: a dataframe with the largest component of each network, the replica and the step
    """

    components = []
    for replica in range(len(nets)):
        for step in range(len(nets[0])):
            net = nets[replica][step]
            largest_component = max(nx.connected_components(net), key=len)
            for node in largest_component:
                components.append({'replica': replica, 'step': step, 'node_number': node})

    return pd.DataFrame(components)



def calculate_indicators(data):
    """
    Calculate the standard deviation, skewness, and kurtosis of the heterozygosity data.
    :param data: DataFrame with 'replica', 'step', and 'het' columns.
    :return: DataFrame with the calculated indicators for each step and replica.
    """

    grouped = data.groupby(['replica', 'step'])['het']
    indicators = grouped.agg(['std', 'skew']).reset_index()
    indicators['kurt'] = grouped.apply(pd.Series.kurtosis).values
    indicators.to_csv(f'indicators.csv', index=False)

    return indicators


def plot_het_indicator(cor: pd.DataFrame,
                       indicators: pd.DataFrame,
                       indicator: str,
                       n_samples: int = 10) -> None:
    """
    Creates a combined plot with two y-axes:
      - Left y-axis: overall mean and standard deviation (SD) of heterozygosity (het) across all replicates,
        plus n_samples random individual replica curves for het.
      - Right y-axis: overall mean and SD of the specified indicator across all replicates,
        plus n_samples random individual replica curves for that indicator.

    The x-axis shows the 'step' converted to a percentage (% fragmentation).

    Parameters:
    -----------
    cor : pd.DataFrame
        DataFrame containing columns: ['step', 'het', 'replica', 'node_number'].
    indicators : pd.DataFrame
        DataFrame containing columns: ['replica', 'step', 'std', 'skew', 'kurt'].
    indicator : str
        The indicator to plot (e.g. 'skew', 'std', or 'kurt').
    n_samples : int, optional
        Number of random individual replica curves to plot for each metric. Default is 10.

    Returns:
    --------
    None
        Displays the combined plot.
    """
    # Use a common scaling for the x-axis: convert 'step' to % fragmentation.
    global_max = max(cor['step'].max(), indicators['step'].max())

    # ---------------------------
    # Overall Heterozygosity Stats
    # ---------------------------
    stats_het = cor.groupby('step')['het'].agg(mean='mean', std='std').reset_index()
    stats_het['step_pct'] = stats_het['step'] / global_max * 100

    # ---------------------------
    # Overall Indicator Stats
    # ---------------------------
    stats_ind = indicators.groupby('step')[indicator].agg(mean='mean', std='std').reset_index()
    stats_ind['step_pct'] = stats_ind['step'] / global_max * 100

    # ---------------------------
    # Set up the figure with two y-axes
    # ---------------------------
    plt.style.use('seaborn-v0_8-whitegrid')
    fig, ax1 = plt.subplots(figsize=(12, 8))
    ax2 = ax1.twinx()
    color_het = 'darkorange'
    color_ind = 'forestgreen'

    # Plot overall het mean and SD on left axis (ax1)
    ax1.plot(stats_het['step_pct'], stats_het['mean'],
             color=color_het)
    ax1.fill_between(stats_het['step_pct'],
                     stats_het['mean'] - stats_het['std'],
                     stats_het['mean'] + stats_het['std'],
                     color=color_het, alpha=0.2)

    # Plot overall indicator mean and SD on right axis (ax2)
    ax2.plot(stats_ind['step_pct'], stats_ind['mean'],
             color=color_ind)
    ax2.fill_between(stats_ind['step_pct'],
                     stats_ind['mean'] - stats_ind['std'],
                     stats_ind['mean'] + stats_ind['std'],
                     color=color_ind, alpha=0.2)

    # ---------------------------
    # Plot individual sample curves for both het and indicator
    # ---------------------------
    # Identify replicas present in both DataFrames.
    common_replicas = list(set(cor['replica'].unique()).intersection(set(indicators['replica'].unique())))
    selected_replicas = random.sample(common_replicas, n_samples)

    for i, replica in enumerate(selected_replicas):
        # --- Heterozygosity for this replica ---
        rep_cor = cor[cor['replica'] == replica]
        rep_het = rep_cor.groupby('step', as_index=False)['het'].mean()
        rep_het['step_pct'] = rep_het['step'] / global_max * 100

        # --- Indicator for this replica ---
        rep_ind = indicators[indicators['replica'] == replica].sort_values('step')
        rep_ind['step_pct'] = rep_ind['step'] / global_max * 100

        # Plot sample curves with lower opacity.
        ax1.plot(rep_het['step_pct'], rep_het['het'], color=color_het, alpha=0.5)
        ax2.plot(rep_ind['step_pct'], rep_ind[indicator], color=color_ind, alpha=0.5)

    ax1.set_xlabel('% fragmentation', fontsize=36)
    ax1.set_ylabel('Heterozygosity', color=color_het, fontsize=36)
    ax2.set_ylabel("Kurtosis", color=color_ind, fontsize=36)

    ax1.tick_params(axis='y', labelsize=32, labelcolor=color_het)
    ax2.tick_params(axis='y', labelsize=32, labelcolor=color_ind)
    ax1.tick_params(axis='x', labelsize=32)

    plt.savefig(f'./figs/het_{indicator}_signelpop.svg', format='svg')
    plt.show()
