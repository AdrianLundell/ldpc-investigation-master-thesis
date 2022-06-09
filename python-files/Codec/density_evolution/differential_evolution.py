# %% Imports and variables
from time import process_time
import numpy as np
from scipy.optimize import linprog
from scipy.linalg import null_space
import random

R = 0.5
dv = 10
dc = 10

# %% Start point
# Assume degree distributions starting from 2
k = 1/np.arange(2, dc+1)
i = 1/np.arange(2, dv+1)

# x = [rho;lam]

c1 = np.concatenate((k, -(1-R)*i))
c2 = np.concatenate((np.ones(dc-1), np.zeros(dv-1)))
c3 = np.concatenate((np.zeros(dc-1), np.ones(dv-1)))
C = np.stack((c1, c2, c3))

d = np.array([[0], [1], [1]])
lb = np.zeros((dv+dc-2, 1))
f = np.ones((dv+dc-2, 1))

sol = linprog(f, A_ub=None, b_ub=None, A_eq=C, b_eq=d, method='interior-point')
x_0 = sol.x

# %% Check that the initial condition

rho = x_0[:dc-1]
lam = x_0[dc-1:]
R_check = 1-np.matmul(k, rho)/np.matmul(i, lam)
print(np.sum(rho), np.sum(lam))

# %% Find the null space to C.

C_c = null_space(C)

# %% Add degree 1
x_0 = np.concatenate(([0], x_0[:dc-1], [0], x_0[dc-1:]))
D = np.size(C_c, 1)
z_vec = np.zeros((1, D))
C_c = np.vstack((z_vec, C_c[:dc-1, :], z_vec, C_c[dc-1:, :]))

# %% Help functions


def compute_x(x_0, C_c, eta):
    return x_0 + np.matmul(C_c, eta)


def poly_eval(p, x):
    exp = np.arange(0, p.size)
    return p*x**exp


def compute_cost(x, dc):
    cost = 0
    for xi in x:
        if xi < 0:
            cost += 100 - 10*xi

    if cost == 0:
        rho = x[:dc]
        lam = x[dc:]

    return cost

# %% Optimization


Np = 100
gens = 10000
F = 0.5  # Mutation variable. Usually in interval [0.1,1]
Cr = 0.7  # Recombination probability
D = np.size(C_c, 1)


t1_start = process_time()
# x = x_0 + eta*C_c
# Initialize population
eta_Np = np.random.rand(D, Np)
cost_Np = np.zeros(Np)
cost_gen = np.zeros(gens)

# Optimize population over generations
for g in range(gens):
    u_Np = np.copy(eta_Np)
    for i in range(Np):
        # Generate selection parameters according to algorithm
        r0 = random.randint(0, Np-1)
        r1 = random.randint(0, Np-1)
        r2 = random.randint(0, Np-1)
        jrand = random.randint(0, D-1)
        while (r0 == i):
            r0 = random.randint(0, Np-1)
        while (r1 == i or r1 == r0):
            r1 = random.randint(0, Np-1)
        while (r2 == i or r2 == r0 or r2 == r1):
            r2 = random.randint(0, Np-1)

        # Mutation and cross over
        cr_arr = np.random.rand(D)
        for j in range(D):
            if (cr_arr[j] < Cr or j == jrand):
                u_Np[j, i] = eta_Np[j, r0] + F*(eta_Np[j, r1]-eta_Np[j, r2])

        # Selection
        if g == 0:
            x = compute_x(x_0, C_c, eta_Np[:, i])
            cost_Np[i] = compute_cost(x, dc)

        x = compute_x(x_0, C_c, u_Np[:, i])
        cost = compute_cost(x, dc)

        if cost < cost_Np[i]:
            eta_Np[:, i] = u_Np[:, i]
            cost_Np[i] = cost

    cost_gen[g] = np.sum(cost_Np)/Np

t1_stop = process_time()
print("Elapsed time:",
      t1_stop-t1_start)

# %% Get the best solution

idx = np.argmin(cost_Np)
#idx = idx[0]
eta_best = eta_Np[:, idx]
x_best = compute_x(x_0, C_c, eta_best)
rho_best = x_best[:dc]
lam_best = x_best[dc:]

# %%
