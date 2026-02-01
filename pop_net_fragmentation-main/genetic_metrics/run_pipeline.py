import pickle
from funcs_run_pipeline import run_replicates
from networks_generator import make_rgg

######### run pipeline
print("here i start!")
n_nodes = 50  # no. of nodes
n_rep = 100
n_edges = 250
net = "RGG"

# # # create list off nets
nets =  make_rgg(n_graphs=n_rep, n_nodes=n_nodes, target_edges=n_edges,asymmetric=True)
print("nets created")


# run the pipeline for all fragmentation types
rand = run_replicates(nets=nets, frag_key='rand',n_workers=20)
print("1")
pickle_filename = f'{net}, rand_asymmetric.pickle'
with open(pickle_filename, 'wb') as file:
    pickle.dump(rand, file)
del rand

cor = run_replicates(nets=nets, frag_key='cor',n_workers=20)
print("2")
with open(f'{net}, cor_asymmetric.pickle', 'wb') as file:
    pickle.dump(cor, file)
del cor

intr = run_replicates(nets=nets, frag_key='intr',n_workers=20)
print("3")
with open(f'{net}, intr_asymmetric.pickle', 'wb') as file:
    pickle.dump(intr, file)
del intr

reg = run_replicates(nets=nets, frag_key='reg',n_workers=20)
print("4")
with open(f'{net}, reg_asymmetric.pickle', 'wb') as file:
    pickle.dump(reg, file)
del reg

div = run_replicates(nets=nets, frag_key='div',n_workers=20)
print("5")
with open(f'{net}, div_asymmetric.pickle', 'wb') as file:
    pickle.dump(div, file)
del div
#
dist = run_replicates(nets=nets, frag_key='dist',n_workers=20)
print("6")
with open(f'{net}, dist_asymmetric.pickle', 'wb') as file:
    pickle.dump(dist, file)
del dist
#
opt = run_replicates(nets=nets, frag_key='opt',n_workers=20)
print("7")
with open(f'{net}, opt_asymmetric.pickle', 'wb') as file:
    pickle.dump(opt, file)
del opt
#
wrst = run_replicates(nets=nets, frag_key='wrst',n_workers=20)
print("8")
with open(f'{net}, wrst_asymmetric.pickle', 'wb') as file:
    pickle.dump(wrst, file)
del wrst
print("finish")
########## finish pipeline