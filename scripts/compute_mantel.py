### calculate mantel for all for fragmentation types
import pandas as pd

from funcs import load_data
from mantel_plot import plot_mantel, plot_fst_distance
from mantel_test import calculate_mantel_all

### compute mantel correlation and write to csv
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
data = load_data(fragmentation_types)
perms = 999
cor_data = calculate_mantel_all(data, perms, dist_type='euclidean')
### plot final figure for fst-distance correlation
data= pd.read_csv('./csv_new/fst_euclidean_corrrlation.csv')
plot_mantel(data)

### plot single step correlation fst-distance
fragmentation_type = ['wrst']
data = load_data(fragmentation_types)
plot_fst_distance(data, fragmentation_type, steps=[0, 75, 150], output_path='./figs/cor_fst_euclidean.svg')

