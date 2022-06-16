"""
Contains graph algorithms for use in qc-ldpc optimisation
"""
import numpy as np 
import galois
from itertools import combinations
import matplotlib.pyplot as plt

def rk_edge_local_girth_layer(G, current_vn_index, rk, t, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, gcd = False, vn_distances = None, cn_distances = None):
    """
    DFS calculation of the rk-edge local girth based on Algorithm 2 in 
    https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
    """

    for i in range(int(enumerated_cn_max[t])):
        current_cn_index = i  #Reduncy to allow differentiating between i and node_i
        
        if not G.has_cyclical_edge_set((current_cn_index, current_vn_index)):
            enumerated_cn_indexes[t] = current_cn_index
            enumerated_cn_max[t+1] = i

            G.add_cyclical_edge_set(current_cn_index, current_vn_index) 
        
            if gcd:
                new_girth = shortest_cycles_gcd(G, current_cn_index, current_vn_index, vn_distances, cn_distances)
            else:
                new_girth = shortest_cycles(G, current_cn_index, current_vn_index)

            girths[t+1] = min(girths[t], new_girth)

            if max_girth[0] <= girths[t+1]:
                if t == rk-1: #Iterate over 0...r_k-1 rather than 1...rk
                    cn_girths.flat[enumerated_cn_indexes[0:t+1]] = girths[t+1]
                    max_girth[0] = girths[t+1]

                    if new_girth == np.inf:
                        G.remove_cyclical_edge_set(current_cn_index, current_vn_index)
                        break

                else: 
                    if gcd:
                        vn_indexes = list(range(G.n_cn, G.n_nodes))
                        vn_distances = shortest_distances(G, current_vn_index, vn_indexes)
                        
                        cn_distances = {}
                        T = np.arange(G.N)
                        for ci in range(0, G.n_cn, G.N):
                            cn_indexes = list(G.proto_index(ci)*G.N + G.proto_value(ci + T))
                            cn_distances[ci] = shortest_distances(G, ci, cn_indexes)
                    else:
                        vn_distances = None 
                        cn_distances = None
                
                    rk_edge_local_girth_layer(G, current_vn_index, rk, t+1, enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, vn_distances, cn_distances)
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

    #Calculate local girths from all variable nodes and from the circulant for i = 0, N, m-N
    if gcd:
        vn_indexes = list(range(G.n_cn, G.n_nodes))
        vn_distances = shortest_distances(G, current_vn_index, vn_indexes)
        
        cn_distances = {}
        T = np.arange(G.N)
        for ci in range(0, G.n_cn, G.N):
            cn_indexes = list(G.proto_index(ci)*G.N + G.proto_value(ci + T))
            cn_distances[ci] = shortest_distances(G, ci, cn_indexes)
    else:
        vn_distances = None 
        cn_distances = None

    rk_edge_local_girth_layer(G, current_vn_index, rk, t, 
                        enumerated_cn_indexes, enumerated_cn_max, girths, max_girth, cn_girths, gcd, vn_distances, cn_distances)

    return max_girth, cn_girths

def shortest_distances(G, node,  stop_nodes = []):
    """
    Shortest path to all listed stop nodes starting from node. Implemented with a BFS algorithm.
    """
    Q = [node]
    explored = set(Q)
    distance = 0
    distances = {}

    while Q and stop_nodes:
        distance += 1
        adjecent_nodes = []

        for node in Q:
            for adjecent_node in G.get_adjecent(node):
                if not adjecent_node in explored:
                    explored.add(adjecent_node)
                    adjecent_nodes.append(adjecent_node)
                    
                    distances[adjecent_node] = distance

                    if adjecent_node in stop_nodes:
                        distances[adjecent_node] = distance
                        stop_nodes.remove(adjecent_node)

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

def shortest_cycles_gcd(G, cn, vn, vn_distances, cn_distances):
    """Approximation of local edge girth for short cycles"""    
    T = np.arange(G.N)
    cn_indexes = list(G.proto_index(cn)*G.N + G.proto_value(cn + T))
    vn_indexes = list(G.proto_index(vn)*G.N + G.proto_value(vn + T) + G.n_cn)
    distances = shortest_distances(G, vn, cn_indexes.copy())
    
    delta = [distances.get(key, np.inf) + 1 for key in cn_indexes]
    result = delta[0]

    for t in range(1, G.N-1):
        t1 = np.mod(G.N-t, G.N)
        cond1 = delta[t] + delta[t1]
        
        t1 = np.mod(t-cn, G.N)
        t2 = np.mod(-cn, G.N)
        temp = cn_distances.get(cn_indexes[t1], np.inf)
        if temp == np.inf:
            cond2 = np.inf
        else:   
            cond2 = vn_distances.get(vn_indexes[t], np.inf) + temp.get(vn_indexes[t2], np.inf) + 2
            
        cond3 = delta[t]*G.N / np.gcd(G.N, t)

        result = min([cond1, cond2, cond3, result])

    return result

def graph_stats(G):
        """Plot H and girth distribution"""
        H = G.get_H()

        girths = []
        for vn in range(G.n_cn, G.n_nodes):
            cycles = shortest_cycles(G, vn)
            if cycles:
                girths.append(min(cycles.values()))
            else:
                girths.append(-10)
                
        parity_eqs = np.zeros((G.m*G.N,G.n*G.N))
        parity_eqs[:,-G.m*G.N:] = 0.2

        plt.subplot(2,1,1)
        plt.imshow(1 - np.maximum(H, parity_eqs), cmap="gray")
        plt.subplot(2,1,2)
        plt.hist(girths)
        plt.show()
        

def make_invertable(G):
    """
    Reorder the QC_tanner_graph G such that the last m columns of parity equations in the protograph form a matrix invertible in GF(N+1) 
    <==> last n_cn equations in H invertible in GF(2)
    Solution found through exhausive search.
    """
    H = galois.GF2(G.get_H().astype(int))
    success = False
    H2 = galois.GF2(np.zeros((G.n_cn, G.n_cn)).astype(int))

    for i, column_indexes in enumerate(combinations((range(G.n)), G.m)):
        for i, j in enumerate(column_indexes):
            H2[:,i*G.N:i*G.N+G.N] = H[:, j*G.N:j*G.N+G.N]
        
        if not np.linalg.det(H2) == 0:
            success = True
            break

    if not success:
        print("Invertion not possible")
    
    reminding_columns = [i for i in range(G.n) if i not in column_indexes]
    new_order = reminding_columns + list(column_indexes)
    G_invertible = G.reordered(new_order)

    return G_invertible

def to_degree_distribution(vn_polynomial, n_vn):
    assert sum(vn_polynomial) == 1
    D = np.zeros(n_vn)
    x = 0
    for degree, coeff in enumerate(vn_polynomial):
        if coeff:   
            n_vns = int(coeff * n_vn)
            D[x:x+n_vns] = degree + 1 
            x += n_vns 

    return D

def vn_polynomial_repr(vn_polynomial):
    s = ""
    for degree, coeff in enumerate(vn_polynomial):
        if coeff:
            s += f"{coeff}*x^{degree} + "

    return s[:-2]

def strategy1(max_girth, cn_girths, G, vn_index):
    # 0)
    for cn in range(G.n_cn):
        if G.has_cyclical_edge_set((cn, vn_index)):
            cn_girths[cn] == -max_girth
    # 1)
    survivors = np.argwhere(cn_girths == max_girth)
    if survivors.size == 1:
        return int(survivors)

    # 2)
    cn_local_girths = np.zeros(survivors.size)
    for i in range(survivors.size):
        cn_index = int(survivors[i])
        if G.add_cyclical_edge_set(cn_index, vn_index):
            cn_local_girths[i] = shortest_cycles(G, cn_index, vn_index)
            G.remove_cyclical_edge_set(cn_index, vn_index)
    survivors = survivors[cn_local_girths == np.max(cn_local_girths)]
   
    if survivors.size == 1:
        return int(survivors)   
    # 3)
    check_degrees = np.take(G.get_check_degrees(), survivors)
    survivors = survivors[check_degrees == np.min(check_degrees)]
    # 4)
    result = survivors[np.random.randint(survivors.size)]

    return int(result)