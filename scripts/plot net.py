import os

from funcs import load_data
from giant_comp import plot_het_component

from variance import process_variance, plot_variance

directory = "/home/lab2/PycharmProjects/fragmentation"
os.chdir(directory)

#### plot heterozygosity-giant compoent
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
data = load_data(fragmentation_types)
plot_het_component(data, frag_types=fragmentation_types, n_bins=20, output='./het_component.svg')


# ############ calculate and plot variance across nodes in the network
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt','wrst']
data = load_data(fragmentation_types)
#for single frag type use fragmentation_types[x]
var = process_variance(data, fragmentation_types)
plot_variance(var)