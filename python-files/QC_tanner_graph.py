import matplotlib.pyplot as plt 
import numpy as np

class QC_tanner_graph:
    """
    Sparse implementation of a tanner graph with protograph dimensions m x n 
    and a scaling factor N
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
            assert self.n_cn <= j < self.n_nodes, "Variable node index out of bounds."
            self.nodes[i].add(j) 
            self.nodes[j].add(i)

    def remove_edges(self, edges):
        """Remove edges defined as pairs of check nodes and variable nodes from the graph"""
        for i, j in edges:
            assert 0 <= i < self.n_cn, "Check node index out of bounds."
            assert self.n_cn <= j < self.n_nodes, "Variable node index out of bounds."
            self.nodes[i].remove(j) 
            self.nodes[j].remove(i)       

    def get_adjecent(self, node):
        """Return all adjecent nodes"""
        return self.nodes[node]

    def has_edge(self, edge):
        """Returns true if the graph contains the edge (ci, vi)"""
        return edge[1] in self.nodes[edge[0]]

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
        assert self.n_cn <= vn_index < self.n_nodes, "Variable node index out of bounds."
        t = np.arange(self.N)
        check_nodes = np.floor(cn_index/self.N) * self.N + np.mod(cn_index + t, self.N)
        variable_nodes = np.floor((vn_index - self.n_cn)/self.N) * self.N + np.mod(vn_index - self.n_cn + t, self.N) + self.n_cn
      
        self.add_edges(np.stack((check_nodes.astype(int), variable_nodes.astype(int)), axis=-1))

    def remove_cyclical_edge_set(self, cn_index, vn_index):
        """Removes a cyclical edge set pi(ci, vi, N) from the graph"""
        assert 0 <= cn_index < self.n_cn, "Check node index out of bounds."
        assert self.n_cn <= vn_index < self.n_nodes, "Variable node index out of bounds."
        t = np.arange(self.N)
        check_nodes = np.floor(cn_index/self.N) * self.N + np.mod(cn_index + t, self.N)
        variable_nodes = np.floor((vn_index - self.n_cn)/self.N) * self.N + np.mod(vn_index - self.n_cn + t, self.N) + self.n_cn
      
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


def run_tests():
    
    G = QC_tanner_graph(3,5,1)
    G.add_edges(zip([0,0],[3,4]))
    G.plot()

if __name__ == "__main__":
    run_tests()
