"""This script implements the MM-QC-PEGA Algorithm
Source: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708 [algorithm 1]
"""

#%% Imports
import numpy as np 
import matplotlib.pyplot as plt 
import qc_graph_tools as qc
import copy 

# %% BFS-Calculation of the rk-Edge local girth
def rk_edge_local_girth_layer(Gt1, current_vn_index, rk, t, enumerated_cn_indexes, enumerated_cn_max, girths, vn_girth, cn_girths):
    """
    DFS calculation of the rk-edge local girth based on Algorithm 2 in 
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
    """

    for i in range(int(enumerated_cn_max[t])):
        current_cn_index = i  #Reduncy to allow differentiating between i and node_i
        
        if not G.has_edge((current_cn_index, current_vn_index)):
            enumerated_cn_indexes[t] = current_cn_index
            enumerated_cn_max[t+1] = i

            Gt2 = copy.deepcopy(Gt1)
            Gt2.add_cyclical_edge_set(current_cn_index, current_vn_index) 
            girths[t+1] = min(girths[t], qc.local_girth_cn(Gt2, current_cn_index, current_vn_index))

            if vn_girth[0] <= girths[t+1]:
                if t == rk-1: #Iterate over 0...r_k-1 rather than 1...rk
                    cn_girths.flat[enumerated_cn_indexes[0:t+1]] = girths[t+1]
                    vn_girth[0] = girths[t+1]
                else: 
                    rk_edge_local_girth_layer(Gt2, current_vn_index, rk, t+1, enumerated_cn_indexes, enumerated_cn_max, girths, vn_girth, cn_girths)
            else:
                pass        
    
def rk_edge_local_girth(G0, current_vn_index, rk):
    t = 0
    enumerated_cn_indexes = np.zeros(rk+1, dtype=int) #s in article
    enumerated_cn_max = np.zeros(rk+1, dtype=int) #u in article
    girths = np.zeros(rk+1) 
    vn_girth = np.array([-np.inf])
    cn_girths = np.full(G0.n_cn, -np.inf)

    enumerated_cn_max[0] = G0.n_cn
    girths[0] = np.inf 

    rk_edge_local_girth_layer(G0, current_vn_index, rk, t, 
                        enumerated_cn_indexes, enumerated_cn_max, girths, vn_girth, cn_girths)

    return vn_girth, cn_girths


#%%
def strategy1(vn_girth, cn_girths, G, vn_index):
    # 1)
    survivors = np.argwhere(cn_girths == vn_girth)
    # 2)
    cn_local_girths = np.zeros(survivors.size)
    for i in range(survivors.size):
        cn_index = int(survivors[i])
        G.add_cyclical_edge_set(cn_index, vn_index)
        cn_local_girths[i] = qc.local_girth_cn(G, cn_index, vn_index)
        G.remove_cyclical_edge_set(cn_index, vn_index)
    survivors = survivors[cn_local_girths == np.max(cn_local_girths)]
    # 3)
    check_degrees = np.take(G.get_check_degrees(), survivors)
    survivors = survivors[check_degrees == np.min(check_degrees)]
    # 4)
    result = survivors[np.random.randint(survivors.size)]

    return int(result)

#%%
m = 2
n = 3
N = 5
d = 2
r = 1
G = qc.QC_tanner_graph(m, n, N)
G.add_cyclical_edge_set(0,m*N)

for j in range(0,n*N,N):
    current_vn_index = j
    for k in range(1, d+1):
        rk = min(r, d - k +1)

        G_temp = copy.deepcopy(G)
        vn_girth, cn_girths = rk_edge_local_girth(G_temp, current_vn_index, rk)
    
        ci = strategy1(vn_girth, cn_girths, G, current_vn_index)
        G.add_cyclical_edge_set(ci, current_vn_index)

plt.spy(G.get_H())