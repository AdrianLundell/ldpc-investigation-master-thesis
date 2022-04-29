import matplotlib.pyplot as plt 
import numpy as np

class QC_tanner_graph:
    """
    Sparse implementation of a tanner graph with protograph dimensions m x n 
    and a scaling factor N
    
    Check nodes are stored in self.nodes[0:n_cn] and variable nodes are stored
    in self.nodes[n_cn:n_nodes], but external methods should not care about this
    and variable nodes should be indexed from 0,...,n_vn-1
    """

    def __init__(self, m, n, N):
        """Creates a graph with m*N check nodes, n*N variable nodes and no edges"""
        self.n_nodes = (m+n)*N
        self.n_cn = m*N
        self.n_vn = n*N 
        self.N = N  

        self.nodes = [set() for i in range(self.n_nodes)]

    def __repr__(self):
        return f"{self.n_cn} CNs, {self.n_vn} VNs, {sum([len(node) for node in self.nodes])//2} edges"

    def add_edges(self, edges):
        """Add edges defined as pairs of check nodes and variable nodes to the graph"""
        for i, j in edges:
            assert 0 <= i < self.n_cn, "Check node index out of bounds."
            assert 0 <= j < self.n_vn, "Variable node index out of bounds."
            j = j + self.n_cn
            self.nodes[i].add(j) 
            self.nodes[j].add(i)

    def remove_edges(self, edges):
        """Remove edges defined as pairs of check nodes and variable nodes from the graph"""
        for i, j in edges:
            assert 0 <= i < self.n_cn, "Check node index out of bounds."
            assert 0 <= j < self.n_vn, "Variable node index out of bounds."
            j = j + self.n_cn
            self.nodes[i].remove(j) 
            self.nodes[j].remove(i)       

    def get_adjecent(self, node):
        """Returns all adjecent nodes of node index node"""
        return self.nodes[node]

    def get_adjecent_cn(self, cn):
        """Returns adjecent nodes of check node index cn"""
        assert 0 <= cn < self.n_cn, "Check node index out of bounds."
        return self.get_adjecent(cn)

    def get_adjecent_vn(self, vn):
        """Returns adjecent nodes of variable node index vn"""
        assert 0 <= vn < self.n_vn, "Variable node index out of bounds."
        return self.get_adjecent(vn + self.cn)

    def has_edge(self, edge):
        """Returns true if the graph contains the edge (ci, vi)"""
        assert 0 <= edge[0] < self.n_cn, f"Edge non existent, check node index {edge[0]} out of bounds."
        assert 0 <= edge[1] < self.n_nodes, f"Edge non existens ariable node index {edge[1]} out of bounds."
        return edge[1] + self.n_cn in self.nodes[edge[0]]

    def get_check_degrees(self) -> list:
        """Returns the degree of all checknodes of the graph"""
        return [len(self.nodes[i]) for i in range(self.n_cn)]

    def get_H(self):
        """Generates a dense representation of the graph"""
        H = np.zeros((self.n_cn, self.n_vn))
        
        for i, nodes in enumerate(self.nodes[0:self.n_cn]):
            for j in nodes:
                H[i,j-self.n_cn] = 1

        return H
    
    def add_cyclical_edge_set(self, cn_index, vn_index):
        """Adds a cyclical edge set pi(ci, vi, N) to the graph"""
        assert 0 <= cn_index < self.n_cn, "Check node index out of bounds."
        assert 0 <= vn_index < self.n_vn, "Variable node index out of bounds."
        t = np.arange(self.N)
        check_nodes = np.floor(cn_index/self.N) * self.N + np.mod(cn_index + t, self.N)
        variable_nodes = np.floor((vn_index)/self.N) * self.N + np.mod(vn_index + t, self.N)
      
        self.add_edges(np.stack((check_nodes.astype(int), variable_nodes.astype(int)), axis=-1))

    def remove_cyclical_edge_set(self, cn_index, vn_index):
        """Removes a cyclical edge set pi(ci, vi, N) from the graph"""
        assert 0 <= cn_index < self.n_cn, "Check node index out of bounds."
        assert 0 <= vn_index < self.n_vn, "Variable node index out of bounds."
        t = np.arange(self.N)
        check_nodes = np.floor(cn_index/self.N) * self.N + np.mod(cn_index + t, self.N)
        variable_nodes = np.floor(vn_index/self.N) * self.N + np.mod(vn_index + t, self.N)
      
        self.remove_edges(np.stack((check_nodes.astype(int), variable_nodes.astype(int)), axis=-1))

    def plot(self):
        
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

def shortest_path(G, vn_index, cn_index = None):
    """Shortest path length from start to stop"""
    if not cn_index is None:
        assert 0 <= cn_index < G.n_cn, "Check node index out of bounds."
    assert 0 <= vn_index < G.n_vn, "Variable node index out of bounds."

    start_index = vn_index + G.n_cn
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

def run_tests():
    m = 4
    n = 4
    N = 1
    G = QC_tanner_graph(m, n, N)
    a = [0, 0, 0, 1, 1, 2, 2, 3, 3]
    b = [0, 1, 2, 0, 1, 2, 3, 2, 3]
    G.add_edges(zip(b,a))

    print("Minimal cycle from vn1 = ", local_girth_vn(G, 0))
    print("Minimal cycle through cn3, vn1 = ",local_girth_cn(G, 2, 0))
    G.plot()

    G = QC_tanner_graph(m, n, 3)
    G.add_cyclical_edge_set(0, 0)
    G.add_cyclical_edge_set(3, 5)
    G.add_cyclical_edge_set(6, 9)

    G.plot()

if __name__ == "__main__":
    run_tests()
