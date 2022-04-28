"""This script implements the MM-QC-PEGA Algorithm
Source: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708 [algorithm 1]
"""

#%% Imports
import numpy as np 
import matplotlib.pyplot as plt 
from QC_tanner_graph import QC_tanner_graph
import copy 

# %% BFS-Calculation of the rk-Edge local girth
def shortest_path(G, start_index, cn_index = None):
    """Shortest path length from start to stop"""
    Q = [start_index]
    explored = {start_index : []}
    
    while Q:
        adjecent_nodes = []

        for node in Q:
            for adjecent_node in G.get_adjecent(node):
                if not adjecent_node in explored:
                    explored[adjecent_node] = explored[node] + [adjecent_node]
                    adjecent_nodes.append(adjecent_node)
                else:
                    overlap = [n1==n2 for n1, n2 in zip(explored[adjecent_node], explored[node])]
                    if not any(overlap) and len(explored[node]) > 1 :
                        if cn_index == explored[adjecent_node][0] or cn_index == explored[node][0] or cn_index is None:
                            return len(explored[node]) + len(explored[adjecent_node]) + 1
    
        Q = adjecent_nodes

    return np.inf

def local_girth_vn(G, vn_index):
    """Minimum cycle length passing through vn"""
    return shortest_path(G, vn_index)

def local_girth_cn(G, cn_index, vn_index):
    """Minimum cycle length passing through the edge (cn, vn) assuming edge between them is known to exist"""
    result = shortest_path(G, vn_index, cn_index)
    
    return result 

def rk_edge_local_girth_layer(G, current_vn_index, rk, t, enumerated_cn_indexes, enumerated_cn_max, girths, vn_girths, cn_girths):
    """
    DFS calculation of the rk-edge local girth based on Algorithm 2 in 
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
    """

    for i in range(int(enumerated_cn_max[t])):
        current_cn_index = i  #Reduncy to allow differentiating between i and node_i
        
        if not G.has_edge((current_cn_index, current_vn_index)):
            enumerated_cn_indexes[t] = current_cn_index
            enumerated_cn_max[t+1] = i

            G_temp = copy.deepcopy(G)
            G_temp.add_cyclical_edge_set(current_cn_index, current_vn_index) 
            girths[t+1] = min(girths[t], local_girth_cn(G_temp, current_cn_index, current_vn_index))

            if vn_girths[t] <= girths[t+1]:
                if t == rk-1: #Iterate over 0...r_k-1 rather than 1...rk
                    vn_girths[t+1] = girths[t+1]
                    cn_girths[t+1, 0:t+1] = girths[t+1]
                else: 
                    rk_edge_local_girth_layer(G_temp, current_vn_index, rk, t+1, enumerated_cn_indexes, enumerated_cn_max, girths, vn_girths, cn_girths)
            else:
                pass        
    

def rk_edge_local_girth(G, current_vn_index, rk):
    t = 0
    enumerated_cn_indexes = np.zeros(rk+1) #s in article
    enumerated_cn_max = np.zeros(rk+1) #u in article
    girths = np.zeros(rk+1) 
    vn_girths = np.zeros(rk+1) #g in article
    cn_girths = np.zeros((rk+1, G.n_cn))

    enumerated_cn_max[0] = G.n_cn
    girths[0] = np.inf 
    vn_girths[0] = -np.inf
    cn_girths[0,:] = -np.inf

    rk_edge_local_girth_layer(G, current_vn_index, rk, t, 
                        enumerated_cn_indexes, enumerated_cn_max, girths, vn_girths, cn_girths)

    return vn_girths[-1], cn_girths[-1, :]

# m = 2
# n = 3
# N = 5
# G = QC_tanner_graph(m, n, N)

# M = 1
# cns = np.random.randint(0, m*N, M)
# vns = np.random.randint(m*N, (m+n)*N, M)

# for cn, vn in zip(cns, vns):
#     G.add_cyclical_edge_set(cn, vn)

# plt.spy(G.get_H())
# print(G)
# plt.show()
# print(rk_edge_local_girth(G, m*N, 2))

# plt.spy(G.get_H())
# plt.show()

m = 4
n = 4
N = 1
G = QC_tanner_graph(m, n, N)
a = [4, 4, 4, 5, 5, 6, 6, 7, 7]
b = [0, 1, 2, 0, 1, 2, 3, 2, 3]
G.add_edges(zip(b,a))

print(local_girth_vn(G, 4))
print(local_girth_cn(G, 2, 4))

m = 4
n = 4
N = 1
G = QC_tanner_graph(m, n, N)
a = [4, 4, 4, 5, 5, 6, 6, 7, 7]
b = [0, 1, 2, 0, 1, 2, 3, 0, 3]
G.add_edges(zip(b,a))

print(local_girth_vn(G, 4))
print(local_girth_cn(G, 2, 4))

#%% r-EDge M-QC-PEGA
def strategy1(metrics):
    """Select the optimal cn according to selection strategy 1"""
    max_d_indexes = np.argwhere(metrics[:,0] == np.max(metrics[:,0]))
    if max_d_indexes.size == 1:
        return int(max_d_indexes)

    max_ace_indexes = np.argwhere(metrics[:,1] == np.max(metrics[max_d_indexes,1]))
    if max_ace_indexes.size == 1:
        return int(max_ace_indexes)

    min_distance_indexes = np.argwhere(metrics[:,2] == np.min(metrics[max_ace_indexes,2]))   
    if min_distance_indexes.size == 1:
        return int(min_distance_indexes)

    random_choice_index = min_distance_indexes[np.random.randint(min_distance_indexes.size)]
    return int(random_choice_index)

#%%
m = 2
n = 3
N = 5
d = 2
r = 1
G = QC_tanner_graph(m, n, N)
G.add_cyclical_edge_set(0,m*N)

def strategy1(vn_girth, cn_girths, G, vn_index):
    # 1)
    survivors = np.argwhere(cn_girths == vn_girth)
    # 2)
    cn_local_girths = np.zeros(survivors.size)
    for i in range(survivors.size):
        cn_index = int(survivors[i])
        G.add_cyclical_edge_set(cn_index, vn_index)
        cn_local_girths[i] = local_girth_cn(G, cn_index, vn_index)
        G.remove_cyclical_edge_set(cn_index, vn_index)
    survivors = survivors[cn_local_girths == np.max(cn_local_girths)]
    # 3)
    check_degrees = np.take(G.get_check_degrees(), survivors)
    survivors = survivors[check_degrees == np.min(check_degrees)]
    # 4)
    result = survivors[np.random.randint(survivors.size)]

    return int(result)


for j in range(0,n*N,N):
    current_vn_index = int(j + m*N)
    for k in range(1, d+1):
        rk = min(r, d - k +1)

        G_temp = copy.deepcopy(G)
        vn_girth, cn_girths = rk_edge_local_girth(G_temp, current_vn_index, rk)
    
        ci = strategy1(vn_girth, cn_girths, G, current_vn_index)
        G.add_cyclical_edge_set(ci, current_vn_index)

plt.spy(G.get_H())
# %%

# %%
