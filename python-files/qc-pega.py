"""This script implements the QC-PEGA Algorithm
Source: https://arxiv.org/pdf/1605.05123.pdf [algorithm 2 with ACE metric]
"""

#%% Prerequsits
from sre_constants import CHCODES
import numpy as np 
import matplotlib.pyplot as plt 

#%%
class Tanner_graph:

    def __init__(self, n_vn, n_cn) -> None:
        
        self.n_nodes = n_cn + n_vn
        self.n_cn = n_cn 
        self.n_vn = n_vn 
        self.nodes = [set() for i in range(self.n_nodes)]

    def __repr__(self):
        return f"Graph with {self.n_cn} CNs, {self.n_vn} VNs and {sum([len(node) for node in self.nodes])//2} edges"

    def add_edges(self, variable_nodes, check_nodes) -> None:

        assert len(variable_nodes) == len(check_nodes), "Variable nodes and check nodes must be of same length"

        for i, j in zip(check_nodes, variable_nodes):
            assert 0 <= i < self.n_cn, "Check node index out of bounds."
            assert 0 <= j < self.n_vn, "Variable node index out of bounds."
            self.nodes[i].add(j+self.n_cn) 
            self.nodes[j+self.n_cn].add(i)

    def get_check_degrees(self) -> list:
        return [len(self.nodes[i]) for i in range(self.n_cn)]

    def compute_metric(self, j):
        """Compute distance and ACE to one variable node for all checknodes"""
        
        assert j < self.n_vn, "Variable node index out of bounds"
        
        not_explored = {node : 0 for node in range(self.n_nodes)}
        distances = [np.inf for  i in range(self.n_nodes)]
        #aces = [np.inf for i in range(self.n_vn) ]
        
        Q = [j + self.n_cn]
        adjecent_nodes = set()        
        while Q:  
            node_index = Q.pop(0)
            adjecent_nodes = adjecent_nodes | (self.nodes[node_index] & set(not_explored))
            distances[node_index] = not_explored.pop(node_index)

            if not Q:
                for key in not_explored:
                    not_explored[key] += 1                    

                Q = list(adjecent_nodes)
                adjecent_nodes = set()

        D = np.array(distances[0:self.n_cn])
        ace = np.zeros(self.n_cn)
        degree = np.array(self.get_check_degrees())

        return np.stack((D, ace, degree))

    def get_H(self):
        H = np.zeros((self.n_cn, self.n_vn))
        
        for i, nodes in enumerate(self.nodes[0:self.n_cn-1]):
            for j in nodes:
                H[i,j-self.n_cn] = 1

        return H
# #Testing
# G = Tanner_graph(n_vn=2, n_cn=3)
# vns = [0,0,0,1,1]
# cns = [0,1,2,0,1]
# G.add_edges(vns, cns)
# print(G.compute_metric(1))
# print([1, 1, 3, 2, 0])

# G = Tanner_graph(n_vn=3, n_cn=3)
# vns = [0,1,1,2,2]
# cns = [0,0,1,1,2]
# G.add_edges(vns, cns)
# print(G.compute_metric(0))
# print([1, 3, 5, 0, 2, 4])


# G = Tanner_graph(n_vn=3, n_cn=3)
# vns = [0]
# cns = [0]
# G.add_edges(vns, cns)
# print(G.compute_metric(0))
# print([1, np.inf, np.inf, 0, np.inf, np.inf])

#%% 
def select_cn(metrics):
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

def compute_edges(i, j, N, n):
    """
    Computes the set of edges to add to the matrix G based on selected check node i for varible node j.
    Returns edges as pairs of node indexes in a E x 2 numpy array for E edges.
    """ 
    t = np.arange(N)
    check_nodes = np.floor(i/N) * N + np.mod(i + t, N) #Indexing for check nodes starts at index N*n of the tanner graph
    variable_nodes = np.floor(j/N) * N + np.mod(j + t, N)

    return variable_nodes.astype(int), check_nodes.astype(int)





#%% Run algorithm
m = 10
n = 90
N = 10
D = 3

n_nodes = (m+n)*N
G = Tanner_graph(n*N, m*N)

for j in np.arange(0, n*N, N):
    metrics = np.full((m*N, 3), np.inf)

    for k in range(D):
        i = select_cn(metrics)
        vns, cns = compute_edges(i, j, N, n) 
        G.add_edges(vns, cns)
        metrics = G.compute_metric(j)

print(G)
plt.spy(G.get_H())
# %%

# %%
