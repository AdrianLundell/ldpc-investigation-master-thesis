"""
Contains graph algorithms for use in qc-ldpc optimisation
"""
import numpy as np 
import galois
from itertools import combinations


def rk_edge_local_girth_layer(G, current_vn_index, rk, t, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, gcd = False):
    """
    DFS calculation of the rk-edge local girth based on Algorithm 2 in 
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
    """

    for i in range(int(enumerated_cn_max[t])):
        current_cn_index = i  #Reduncy to allow differentiating between i and node_i
        
        # print(i)
        if not G.has_cyclical_edge_set((current_cn_index, current_vn_index)):
            enumerated_cn_indexes[t] = current_cn_index
            enumerated_cn_max[t+1] = i

            G.add_cyclical_edge_set(current_cn_index, current_vn_index) 
            girths[t+1] = min(girths[t], shortest_cycles(G, current_cn_index, current_vn_index))

            if max_girth[0] <= girths[t+1]:
                if t == rk-1: #Iterate over 0...r_k-1 rather than 1...rk
                    cn_girths.flat[enumerated_cn_indexes[0:t+1]] = girths[t+1]
                    max_girth[0] = girths[t+1]
                else: 
                    #Calculate  Fv, Fc,c for GCD
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

def shortest_distances(G, node,  stop_node = None):
    """
    Shortest path to all connected nodes starting from node. Implemented with a BFS algorithm.

    If a stop node is provided, the algorithm stops and the current distance is returned.
    """
    Q = [node]
    explored = set(Q)
    distance = 0
    distances = {}

    while Q:
        distance += 1
        adjecent_nodes = []

        for node in Q:
            for adjecent_node in G.get_adjecent(node):
                if not adjecent_node in explored:
                    explored.add(adjecent_node)
                    adjecent_nodes.append(adjecent_node)
                    
                    distances[adjecent_node] = distance

                    if adjecent_node == stop_node:
                        return distance

        Q = adjecent_nodes

    return distances

def shortest_cycles(G, node, stop_node = None):
    """
    Shortest cycles starting in node and passing through each adjecent node.

    Performs a BFS search, saving the path taken up to that path. This path is necessarily minimal.
    If a stop node is provided, the algorithm stops when the cycle containing this node is found and this cycle distance is returned.
    """
    Q = [node]
    explored = {node : []}
    distances = {}

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
                        n1 = explored[node][0]
                        n2 = explored[adjecent_node][0]

                        distance = len(explored[node]) + len(explored[adjecent_node]) + 1
                        distances[n1] = distance
                        distances[n2] = distance

                        if stop_node in explored[node]:
                            return distance

        Q = adjecent_nodes
    if stop_node is None:
        return distances
    else:
        return np.inf

# def shortest_cycles_gcd(G, vn, cn, vn_distances, cn_distances):
#     """Approximation of local edge girth for short cycles"""
#     distances = shortest_distances(G, vn)
#     distances = shortest_distances(G, vn)

#     if not distances:
#         return np.inf
    
#     shifted_cn = G.shift(cn, np.arange(G.N))
#     shifted_vn = G.shift(vn, np.arange(G.N)) + G.n_vn
    
#     delta = [distances.get(key, np.inf) + 1 for key in shifted_cn]
#     result = delta[0]

#     for t in range(1, G.N):
#         condition1 = delta[t] + delta[G.N-(t+1)]
#         condition2 = distances.get(shifted_vn[t], np.inf) +  
#         min_distance(G, G.shift(cn, t-cn), G.shift(cn, -cn))
#         cond3 = delta[t]*G.N / np.gcd(G.N, t)

#         result = min([cond1, cond2, cond3, result])

#     return result

def make_invertable(G):
    """
    Reorder the QC_tanner_graph G such that the last m columns of parity equations in the progograph form a matrix invertible in GF(N+1) 
    <==> last n_cn equations in H invertible in GF(2)
    Solution found through exhausive search.
    """
    proto = G.proto
    n_rows,n_cols = proto.shape
    GF = galois.GF(G.N + 1)

    for column_indexes in combinations((range(n_cols)), n_rows):
        gf_mat = GF(np.take(proto, column_indexes, axis=-1))
        
        if not np.linalg.det(gf_mat) == 0:
            break

    reminding_columns = [i for i in range(n_cols) if i not in column_indexes]
    new_order = reminding_columns + list(column_indexes)
    G_invertible = G.reorder(new_order)

    return G_invertible
