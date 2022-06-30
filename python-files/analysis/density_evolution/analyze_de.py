import numpy as np
import matplotlib.pyplot as plt

######################################
# Discrete analysis
######################################

fname =  "../data/discrete_test2.npz"
data = np.load(fname)
population = data["population"]
fitness = data["fitness"]
last_gen = int(data["generation"][0])+1
best_idx = int(data["best_idx"][0])
lam_node = data["lam_node"]
rho_node = data["rho_node"]

def get_degree_distrubutions(i):
    cn_degrees = np.sum(i, 0)
    vn_degrees = np.sum(i, 1)
    rho_node = np.bincount(vn_degrees)[1:]
    rho_edge = np.arange(1, len(rho_node)+1) * rho_node / \
        np.sum(np.arange(1, len(rho_node)+1) * rho_node)
    lam_node = np.bincount(cn_degrees)[1:]
    lam_edge = np.arange(1, len(lam_node)+1) * lam_node / \
        np.sum(np.arange(1, len(lam_node)+1) * lam_node)
    rho_node = rho_node/np.sum(rho_node)
    lam_node = lam_node/np.sum(lam_node)
    return lam_edge, rho_edge


def convert_edge_to_node(p_edge):
    p_node = p_edge/(np.arange(1, len(p_edge)+1)*np.sum(1/np.arange(1, len(p_edge)+1)*p_edge))
    return p_node

best_idx = np.argmax(fitness[last_gen,:])
lam_edge, rho_edge = get_degree_distrubutions(population[best_idx,:,:])
print(f"Lambda edge degeree distribution: \n {lam_edge}")
print(f"Rho edge degeree distribution: \n {rho_edge}")
print(f"Lambda node degeree distribution: \n {lam_node}")
print(f"Rho node degeree distribution: \n {rho_node}")

p_node = convert_edge_to_node(lam_edge)
print(f"Lambda node degeree distribution: \n {p_node}")

######################################
# Continuous analysis
######################################
fname = "../data/continuous_test2.npz"
data = np.load(fname)
population = data["population"]
fitness = data["fitness"]
last_gen = int(data["generation"][0])+1
best_idx = int(data["best_idx"][0])
lam_node = data["lam_node"]
rho_node = data["rho_node"]
best_rber = data["best_rber"][0]
print(f"Lambda node degeree distribution: \n {lam_node}")
print(f"Rho node degeree distribution: \n {rho_node}")
print(f"Best RBER: {best_rber}")
