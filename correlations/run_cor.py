### calculate mantel for all for fragmentation types
import pandas as pd

from funcs import load_data
from cor_plot import plot_mantel, plot_fst_distance
from cor_test import calculate_qap_all

## compute mantel correlation and write to csv
fragmentation_types = ['rand', 'cor', 'intr', 'reg', 'dist', 'div', 'opt', 'wrst']
# data = load_data(fragmentation_types)
perms = 1000
# cor_data = calculate_qap_all(data, perms, dist_type='resistance')
# plot final figure for fst-distance correlation
data= pd.read_csv('fst_resistance.csv')
plot_mantel(data)

### plot single step correlation fst-distance
# fragmentation_type = ['wrst']
#
# plot function loads data internally; just pass fragmentation_type list
# plot_fst_distance(fragmentation_type, steps=[0,75,150], dist_type='mfpt')

