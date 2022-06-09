"""
This file validates the implementation of the density evolution by replicating the results of the simple bec-channel expression.

Expected 
(3, 6) 0.4294
(4, 8) 0.3834
(5, 10) 0.3415
(6, 12) 0.3074
(7, 14) 0.2798
"""
import time 
import numpy as np 
import density_evolution 
import matplotlib.pyplot as plt

t0 = time.time()
n_grid = 512

ensamble = (3, 6)
lambda_coeffs = [0]*ensamble[0]
lambda_coeffs[-1] = 1
rho_coeffs = [0]*ensamble[1]
rho_coeffs[-1] = 1

min = 0.1
max = 1

f_grid = np.linspace(-30, 30, n_grid, endpoint = False)
f_step = abs(f_grid[1] - f_grid[0])
max_val = -np.log(np.tanh(f_step))
g_grid = np.linspace(0, max_val, n_grid)

def eval(x):
    cdf = np.zeros(n_grid)
    cdf[n_grid//2:] = x
    cdf[-1] = 1
    return density_evolution.symmetric_density_evolution(cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, plot = False)

result = density_evolution.bisection_search(min, max, eval)

print(result)
print(time.time()-t0)
