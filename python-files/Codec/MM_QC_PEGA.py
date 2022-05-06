"""This script implements the MM-QC-PEGA Algorithm
Source: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
"""

#%% Prerequisites
import numpy as np 
import matplotlib.pyplot as plt 
import Tanner_graph
import Graph_algorithms
import galois 

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
m = 2
n = 10
N = 16
r = 1  
G = Tanner_graph.QC_tanner_graph(m, n, N)
D = np.full(G.n_vn,2)
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

parity_eqs = np.zeros((m*N,n*N))
parity_eqs[:,-m*N:] = 0.2

plt.subplot(2,1,1)
plt.imshow(1 - np.maximum(G.get_H(), parity_eqs), cmap="gray")
plt.subplot(2,1,2)
plt.imshow(1 - np.maximum(G_reordered.get_H(), parity_eqs), cmap="gray")
plt.show()

GF = galois.GF(2)
H = GF(G.get_H().astype(int))
H_reordered = GF(G_reordered.get_H().astype(int))

try:
    print(np.linalg.det(H[:,-m*N:]))
    np.linalg.inv(H[:,-m*N:])
except: 
    print("G invertion failed")

try:
    print(np.linalg.det(H_reordered[:,-m*N:]))
    np.linalg.inv(H_reordered[:,-m*N:])
except: 
    print("G reoredered invertion failed")
#G.save("test.qc")
# %%
