"""This script implements the MM-QC-PEGA Algorithm
Source: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=8241708
"""

#%% Prerequisites
import numpy as np 
import matplotlib.pyplot as plt 
import Tanner_graph
import Graph_algorithms
import galois 
import datetime
import time 

def vn_distribution(vn_polynomial, n_vn):
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
            cn_local_girths[i] = Graph_algorithms.shortest_cycles(G, cn_index, vn_index)
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

#%% MM-QC-PEGA algorithm
def mm_qc_pega(r,m,n,N, vn_polynomial_coeffs, gcd = False, seed = None):
    
    header = f"""
===================================================================
Creating ({m},{n},{N}) code with the {r}-edge QC-PEGA algorithm for 
variable node edge distribution {vn_polynomial_repr(vn_polynomial_coeffs)} 
===================================================================
"""
    print(header)
    print("Initialising...")
    
    G = Tanner_graph.QC_tanner_graph(m, n, N)
    D = vn_distribution(vn_polynomial_coeffs, n*N)
    if not seed == None:
        np.random.seed(seed)
    t0 = time.time()

    for current_vn_index in range(0,n*N,N):
        d = D[current_vn_index]
        for k in range(1, int(d+1)):
            rk = int(min(r, d - k +1))
            max_girth, cn_girths = Graph_algorithms.rk_edge_local_girth(G, current_vn_index, rk, gcd = gcd)
        
            ci = strategy1(max_girth, cn_girths, G, current_vn_index)
            G.add_cyclical_edge_set(ci, current_vn_index)

        dt = time.time() - t0
        completed_progress = (current_vn_index+N)/(n*N)
        time_left = dt/completed_progress - dt        
        status = f"{100*completed_progress:.2f}% completed, elapsed time {int(dt//60)}:{dt% 60:.2f}s. Approximate time left: {int(time_left//60)}:{time_left % 60:.2f}"
        print(status, end='\r', flush=True)
        

    print("")
    print(f"Done. Time taken: {int(dt//60)} minutes, {dt % 60:.2f} seconds.")

    return G

#%%
name = "Test"
vn_polynomial = []
r = 2
m = 9
n = 73
N = 10
G = mm_qc_pega(r, m, n, N, vn_polynomial, False, 0)
Graph_algorithms.graph_stats(G)

#%% Ensure invertible party equations
H = galois.GF2(G.get_H().astype(int))
assert np.linalg.matrix_rank(H) == G.n_cn, "Bad matrix rank"
G_reordered = Graph_algorithms.make_invertable(G)

# %% Save file
metadata = f"""
# QC-code generated with the r-edge metric constrained QC-PEG algorithm.
# 
# Name   : {name}
# D      : {vn_polynomial_repr(vn_polynomial)}
# R      : {r}
# Metric : Minimum distance
# Date   : {datetime.datetime.now()}
"""
G_reordered.save(f"{name}.qc")
f = open(f"{name}.qc", 'a')
f.write(metadata)
f.close()