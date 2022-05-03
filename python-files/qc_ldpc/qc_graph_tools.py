import matplotlib.pyplot as plt 
import numpy as np

class QC_tanner_graph:
    """
    Sparse implementation of a tanner graph with protograph dimensions m x n 
    and a scaling factor N
    
    Check nodes are stored in self.nodes[0:n_cn] and variable nodes are stored
    in self.nodes[n_cn:n_nodes], but external methods should not care about this
    and variable nodes should be indexed from 0,..., n_vn-1 externally.
    """

    def __init__(self, m, n, N):
        """Creates a graph with m*N check nodes, n*N variable nodes and no edges"""
        assert m > 0 and n > 0 and N > 0, "m, n and N must be positive integers."
        
        self.n_nodes = int((m+n)*N)
        self.n_cn = int(m*N)
        self.n_vn = int(n*N) 
        self.N = int(N)  

        self.nodes = [set() for i in range(self.n_nodes)]

    def __repr__(self):
        return f"{self.n_cn} CNs, {self.n_vn} VNs, {sum([len(node) for node in self.nodes])//2} edges"

    def assert_edge(self, edge):
        assert 0 <= edge[0] < self.n_cn, f"Edge non existent, check node index {edge[0]} out of bounds."
        assert 0 <= edge[1] < self.n_nodes, f"Edge non existent, variable node index {edge[1]} out of bounds."

    def assert_vn_node(self, vn):
        assert 0 <= vn < self.n_vn, "Variable node index out of bounds."

    def assert_cn_node(self, cn):
        assert 0 <= cn < self.n_cn, "Check node index out of bounds."
    

    def add_edges(self, edges):
        """Add edges defined as pairs of check nodes and variable nodes to the graph"""
        for edge in edges:
            self.assert_edge(edge)
            self.nodes[edge[0]].add(edge[1] + self.n_cn) 
            self.nodes[edge[1] + self.n_cn].add(edge[0])

    def remove_edges(self, edges):
        """Remove edges defined as pairs of check nodes and variable nodes from the graph"""
        for edge in edges:
            assert self.has_edge(edge), "Cannot remove non existent edge"
            self.nodes[edge[0]].remove(edge[1] + self.n_cn) 
            self.nodes[edge[1] + self.n_cn].remove(edge[0])     

    def get_adjecent(self, node):
        """Returns all adjecent nodes of node index node"""
        return self.nodes[int(node)]

    def get_adjecent_cn(self, cn):
        """Returns adjecent nodes of check node index cn"""
        self.assert_cn_node(cn)
        return self.get_adjecent(cn)

    def get_adjecent_vn(self, vn):
        """Returns adjecent nodes of variable node index vn"""
        self.assert_vn_node(vn)
        return self.get_adjecent(vn + self.n_cn)

    def has_edge(self, edge):
        """Returns true if the graph contains the edge (ci, vi)"""
        self.assert_edge(edge)
        return edge[1] + self.n_cn in self.nodes[edge[0]]

    def get_check_degrees(self) -> list:
        """Returns the degree of all check nodes of the graph"""
        return [len(self.nodes[i]) for i in range(self.n_cn)]

    def get_H(self):
        """Generates a dense representation of the graph"""
        H = np.zeros((self.n_cn, self.n_vn))
        
        for i, nodes in enumerate(self.nodes[0:self.n_cn]):
            for j in nodes:
                H[i,j-self.n_cn] = 1

        return H

    def shift(self, node, t):
        """Shifts the node index by t cyclically"""
        return np.floor(node/self.N) * self.N + np.mod(node + t, self.N)

    def add_cyclical_edge_set(self, cn_index, vn_index):
        """Adds a cyclical edge set pi(ci, vi, N) to the graph"""
        self.assert_edge((cn_index, vn_index))
        t = np.arange(self.N)
        check_nodes = self.shift(cn_index, t)
        variable_nodes = self.shift(vn_index, t)
      
        self.add_edges(np.stack((check_nodes.astype(int), variable_nodes.astype(int)), axis=-1))

    def remove_cyclical_edge_set(self, cn_index, vn_index):
        """Removes a cyclical edge set pi(ci, vi, N) from the graph"""
        self.assert_edge((cn_index, vn_index))
        t = np.arange(self.N)
        check_nodes = np.floor(cn_index/self.N) * self.N + np.mod(cn_index + t, self.N)
        variable_nodes = np.floor(vn_index/self.N) * self.N + np.mod(vn_index + t, self.N)
      
        self.remove_edges(np.stack((check_nodes.astype(int), variable_nodes.astype(int)), axis=-1))

    def local_girth(self, edge, gcd = False):
        """Calculates the local girth of edge (cn, vn), accurately or with gcd-approximation"""
        self.assert_edge(edge)
        if gcd:
            return shortest_cycles_gcd(self, edge[1], edge[0]) 
        else:
            return shortest_cycles(self, edge[1], edge[0])


    def plot(self):
        """Graphical representation of the tanner graph"""
        width = 500
        height = 250
        border = 100
        vn_coords = np.stack((np.linspace(0, width, self.n_vn), np.full(self.n_vn, height)))
        cn_coords = np.stack((np.linspace(0, width, self.n_cn), np.full(self.n_cn, 0)))
        
        #Init figure
        plt.figure()
        plt.xlim(-border, width+border)
        plt.ylim(-border, height + border)
        
        #Plot nodes
        plt.scatter(vn_coords[0,:],vn_coords[1,:], s = 40, c = "black", marker = "o")
        plt.scatter(cn_coords[0,:],cn_coords[1,:], s = 40, c = "black", marker = "s")
        
        #Plot edges
        for cn, vns in enumerate(self.nodes[:self.n_cn]):
            for vn in vns:
                vn_x = vn_coords[0,vn-self.n_cn]
                vn_y = vn_coords[1,vn-self.n_cn]
                cn_x = cn_coords[0, cn]    
                cn_y = cn_coords[1, cn]    
                plt.plot([vn_x, cn_x], [vn_y, cn_y], c = "black")

        plt.show()

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
