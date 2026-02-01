



###### make a csv file for singlepop analysis using the earlywarning R package
# (10 nodes with the highest heterozygosity in the last step of each replica)
# read the pickle file with RGG (d-0.6) data
# with open(f'cor_d0.6_r1000.pickle', 'rb') as file:
#     cor = pickle.load(file)
# print('finish')
# frag = 'cor_d0.6_r1000'
# export_het_csv(cor, frag)



###### make data for metapop analysis, get the het for the largest component
with open('./pickles/cor_d0.6_r1000.pickle', 'rb') as file:
    raw = pickle.load(file)

length, nets, het_dist, het_mean, fst_dist = raw
het = assign_node_numbers(het_dist)
components = get_largest_component(nets)
component_data = pd.merge(het, components, on=['replica', 'step', 'node_number'])
component_data = component_data.sort_values(by=['replica', 'step', 'node_number'])
component_data.reset_index().to_csv('cor_d0.6_r1000_component.csv', index=False)


#### calculate indicators for metapop (largest component)
# cor = pd.read_csv('cor_d0.6_r1000_component.csv')
# indicators = calculate_indicators(cor)
# indicators.to_csv('indicators_metapop.csv', index=False)



###### plot het+indicators of metapopilation data-change y label
# random.seed(1)
# cor = pd.read_csv('cor_d0.6_r1000_het.csv')
# indicators = pd.read_csv('indicators_metapop.csv')
# plot_het_indicator(cor, indicators, indicator='kurt', n_samples=10)

###### plot het+indicators of single population data-change y label
# random.seed(1)
# cor = pd.read_csv('cor_d0.6_r1000_het.csv')
# indicators = pd.read_csv('indicators_singlepop_25.csv')
# plot_het_indicator(cor, indicators, indicator='kurt', n_samples=10)