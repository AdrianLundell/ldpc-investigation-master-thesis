"""This script implements the MM-QC-PEGA Algorithm
Source: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
"""

#%% Prerequisites
import numpy as np 
import matplotlib.pyplot as plt 
import Tanner_graph
import Graph_algorithms

def strategy1(max_girth, cn_girths, G, vn_index):
    # 0)
    for cn in range(G.n_cn):
        if G.has_cyclical_edge_set((cn, vn_index)):
            cn_girths[cn] == -max_girth
    # 1)
    survivors = np.argwhere(cn_girths == max_girth)
    # 2)
    cn_local_girths = np.zeros(survivors.size)
    for i in range(survivors.size):
        cn_index = int(survivors[i])
        if G.add_cyclical_edge_set(cn_index, vn_index):
            cn_local_girths[i] = Graph_algorithms.shortest_cycles(G, cn_index, vn_index)
            G.remove_cyclical_edge_set(cn_index, vn_index)
    survivors = survivors[cn_local_girths == np.max(cn_local_girths)]
    # 3)
    check_degrees = np.take(G.get_check_degrees(), survivors)
    survivors = survivors[check_degrees == np.min(check_degrees)]
    # 4)
    result = survivors[np.random.randint(survivors.size)]

    return int(result)

#%% MM-QC-PEGA algorithm
m = 8
n = 16
N = 5
r = 1  
G = Tanner_graph.QC_tanner_graph(m, n, N)
D = np.full(G.n_vn,3)
np.random.seed(0)

for j in range(0,n*N,N):
    current_vn_index = j
    d = D[j]

    for k in range(1, int(d+1)):
        rk = min(r, d - k +1)

        max_girth, cn_girths = Graph_algorithms.rk_edge_local_girth(G, current_vn_index, rk, gcd = False)
    
        ci = strategy1(max_girth, cn_girths, G, current_vn_index)
        G.add_cyclical_edge_set(ci, current_vn_index)

#%% Invertible check
G_reordered = Graph_algorithms.make_invertable(G)

plt.spy(G.get_H())
plt.show()
plt.spy(G_reordered.get_H())
plt.show()

#G.save("test.qc")