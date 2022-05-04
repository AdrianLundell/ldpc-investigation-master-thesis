"""This script implements the MM-QC-PEGA Algorithm
Source: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708 [algorithm 1]
"""

#%% Imports
import numpy as np 
import matplotlib.pyplot as plt 
import qc_graph_tools as qc
import copy 
from itertools import combinations

# %% BFS-Calculation of the rk-Edge local girth
def rk_edge_local_girth_layer(G, current_vn_index, rk, t, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, gcd = False):
    """
    DFS calculation of the rk-edge local girth based on Algorithm 2 in 
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
    """

    for i in range(int(enumerated_cn_max[t])):
        current_cn_index = i  #Reduncy to allow differentiating between i and node_i
        
        # print(i)
        if not G.has_edge((current_cn_index, current_vn_index)):
            enumerated_cn_indexes[t] = current_cn_index
            enumerated_cn_max[t+1] = i

            G.add_cyclical_edge_set(current_cn_index, current_vn_index) 
            girths[t+1] = min(girths[t], qc.shortest_cycles(G, current_cn_index, current_vn_index))

            if max_girth[0] <= girths[t+1]:
                if t == rk-1: #Iterate over 0...r_k-1 rather than 1...rk
                    cn_girths.flat[enumerated_cn_indexes[0:t+1]] = girths[t+1]
                    max_girth[0] = girths[t+1]
                else: 
                    #Calculate  Fv, Fc,c fopr GCD
                    rk_edge_local_girth_layer(G, current_vn_index, rk, t+1, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths)
            else:
                pass        

            G.remove_cyclical_edge_set(current_cn_index, current_vn_index)

def rk_edge_local_girth(G, current_vn_index, rk, gcd = False):
    """
    Calculate the maximum girth possible when adding an edge from current_vn_index to each check node, with a look-ahead depth of rk. 
    """
    t = 0
    enumerated_cn_indexes = np.zeros(rk+1, dtype=int) #s in article
    enumerated_cn_max = np.zeros(rk+1, dtype=int) #u in article
    girths = np.zeros(rk+1) 
    max_girth = np.array([-np.inf])
    cn_girths = np.full(G.n_cn, -np.inf)

    enumerated_cn_max[0] = G.n_cn
    girths[0] = np.inf 

    rk_edge_local_girth_layer(G, current_vn_index, rk, t, 
                        enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, gcd)

    return max_girth, cn_girths


#%%
def strategy1(max_girth, cn_girths, G, vn_index):
    # 1)
    survivors = np.argwhere(cn_girths == max_girth)
    # 2)
    cn_local_girths = np.zeros(survivors.size)
    for i in range(survivors.size):
        cn_index = int(survivors[i])
        G.add_cyclical_edge_set(cn_index, vn_index)
        cn_local_girths[i] = qc.shortest_cycles(G, cn_index, vn_index)
        G.remove_cyclical_edge_set(cn_index, vn_index)
    survivors = survivors[cn_local_girths == np.max(cn_local_girths)]
    # 3)
    check_degrees = np.take(G.get_check_degrees(), survivors)
    survivors = survivors[check_degrees == np.min(check_degrees)]
    # 4)
    result = survivors[np.random.randint(survivors.size)]

    return int(result)

#%% GCD-approximation of local cn_girth


#%%
# D = np.zeros(G.n_vn)
# D[:264] = 2
# D[264:264+192] = 3
# D[264+192:] = 6

m = 8
n = 16
N = 10
r = 1  
G = qc.QC_tanner_graph(m, n, N)
D = np.full(G.n_vn,3)
#For same result each time
np.random.seed(0)

for j in range(0,n*N,N):
    current_vn_index = j
    d = D[j]

    for k in range(1, int(d+1)):
        # print("k: ", k)
        rk = min(r, d - k +1)

        #Calculate Fv, Fc,c for GCD
        max_girth, cn_girths = rk_edge_local_girth(G, current_vn_index, rk, gcd = False)
    
        ci = strategy1(max_girth, cn_girths, G, current_vn_index)
        G.add_cyclical_edge_set(ci, current_vn_index)
        G.add_to_proto(ci, current_vn_index)

#%%
plt.spy(G.get_H())
plt.show()
np.linalg.inv(G.get_H()[:,:m*N])

#%%
A = G.get_H()

n_rows,n_cols = A.shape

ones = np.zeros(n_rows)
free_ones = np.ones(n_rows)
a = np.zeros(n_rows)
columns = []

for i, column_indexes in enumerate(combinations((range(n_cols)), n_rows)):
    det = np.linalg.det(np.take(A, column_indexes, axis=-1))
    #print(det)

    if not det == 0:
        print(column_indexes)
        break


#%%
G.save("test.qc")
G = qc.QC_tanner_graph.read("/home/adrianlundell/ldpc-investigation-master-thesis/python-files/qc_ldpc/test.qc")
plt.spy(G.get_H())
plt.show()
# %%
for i in range(G.n_vn):
    print(qc.shortest_cycles(G, i + G.n_cn))
    print(len(G.nodes[i + G.n_cn]))

#%%
print(G.proto)
# %%
