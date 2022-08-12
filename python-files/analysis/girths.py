#%%
import sys 
import numpy as np 

sys.path.append('../../mm_qc_pega')
import graphs
import peg_utils
# %%

G = graphs.QC_tanner_graph.read("data/optimal_soft.qc")
#peg_utils.graph_stats(G)


# %%
vn_nodes = np.bincount(G.get_var_degrees())[1:]
cn_nodes = np.bincount(G.get_check_degrees())[1:]

vn_nodes = vn_nodes / np.sum(vn_nodes)
cn_nodes = cn_nodes / np.sum(cn_nodes)

vn_edges = vn_nodes * np.arange(len(vn_nodes))
cn_edges = cn_nodes * np.arange(len(cn_nodes))

vn_edges = vn_edges / np.sum(vn_edges)
cn_edges = cn_edges / np.sum(cn_edges)

print(vn_edges, cn_edges)
# %%
