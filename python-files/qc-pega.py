"""This script implements the MM-QC-PEGA Algorithm
Source: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708 [algorithm 1]
"""

#%% Imports
import numpy as np 
import matplotlib.pyplot as plt 
from QC_tanner_graph import QC_tanner_graph

# %% DFS-Calculation of the rk-Edge local girth

def shortest_path(G, start_index, stop_index):
    """Shortest path length from start to stop"""
    Q = [start_index]
    explored = set(Q)
    length = 0

    while Q:
        length += 1
        adjecent_nodes = []

        for node in Q:
            for adjecent_node in G.get_adjecent(node):
                if adjecent_node == stop_index and length > 1:
                    return length 
                if not adjecent_node in explored:
                    explored.add(adjecent_node)
                    adjecent_nodes.append(adjecent_node)
        Q = adjecent_nodes

    raise Exception("Path not found")

def local_girth_vn(G, vn_index):
    """Minimum cycle length passing through vn"""
    return shortest_path(G, vn_index, vn_index)

def local_girth_cn(G, cn_index, vn_index):
    """Minimum cycle length passing through the edge (cn, vn) assuming edge between them is known to exist"""
    G.remove_edges([(cn_index, vn_index)])
    result = shortest_path(G, cn_index, vn_index) + 1
    G.add_edges([(cn_index, vn_index)])
    
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

            G.add_cyclical_edge_set(current_cn_index, current_vn_index) 
            girths[t+1] = min(girths[t], local_girth_cn(G, current_cn_index, current_vn_index))

            if vn_girths[t] <= cn_girths[t+1, current_cn_index]:
                if t == rk-1: #Iterate over 0...r_k-1 rather than 1...rk
                    vn_girths[t+1] = girths[t+1]
                    cn_girths[t+1, t] = girths[t+1]
                else: 
                    rk_edge_local_girth_layer(G, current_vn_index, rk, t+1, enumerated_cn_indexes, enumerated_cn_max, girths, vn_girths, cn_girths)
                    return 

def rk_edge_local_girth(G, current_vn_index, rk):
    t = 0
    enumerated_cn_indexes = np.zeros(rk) #s in article
    enumerated_cn_max = np.zeros(rk) #u in article
    girths = np.zeros(rk)
    vn_girths = np.zeros(rk) #g in article
    cn_girths = np.zeros((rk, G.n_cn))

    enumerated_cn_max[0] = G.n_cn
    girths[0] = np.inf 
    vn_girths[0] = -np.inf
    cn_girths[0,:] = -np.inf

    rk_edge_local_girth_layer(G, current_vn_index, rk, t, 
                        enumerated_cn_indexes, enumerated_cn_max, girths, vn_girths, cn_girths)

    return vn_girths[-1]

m = 4
n = 8
N = 4
G = QC_tanner_graph(m, n, N)

M = 30
cns = np.random.randint(0, m*N, M)
vns = np.random.randint(m*N, (m+n)*N, M)
for cn, vn in zip(cns, vns):
    G.add_cyclical_edge_set(cn, vn)

#plt.spy(G.get_H())
#print(G)

print(rk_edge_local_girth(G, m*N + 2, 3))

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

