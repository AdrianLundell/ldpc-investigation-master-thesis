"""
This file puts together the methods from discretize_sigma, sampled_de and symmetric_de to 
show correctness of the algorithms by comparing sampled and theoretical distributions.
"""
from math import gamma
import sampled_de
import symmetric_de
import discretize_sigma

import numpy as np 
import matplotlib.pyplot as plt

n_grid = 8192
sigma, p0, bins = discretize_sigma.compute_pdf(0.01, 0.5, n_grid, 30)
f_grid, g_grid, pdf = discretize_sigma.create_pdf(p0, bins, n_grid)
cdf = sampled_de.to_cdf(pdf)
n_samples = 10**5//2
samples = sampled_de.sample(cdf, f_grid, n_samples)
lambda_coeffs = np.array([0, 0, 1])
rho_coeffs = np.array([0, 0, 0, 0, 0, 1])

#Compare initial distributions
#sampled_de.plot_samples(samples, f_grid, n_samples)
# plt.plot(f_grid, cdf)
# plt.show()

#Compare gamma
g1 = sampled_de.gamma(samples)
g2 = symmetric_de.gamma(cdf, f_grid, g_grid)
sampled_de.plot_samples(g1[0], g_grid, n_samples)
sampled_de.plot_samples(g1[1], g_grid, n_samples)
plt.plot(g_grid, g2[0,:])
plt.plot(g_grid, g2[1,:])
plt.show()

#Compare rho
r1 = sampled_de.rho(g1, rho_coeffs, n_samples)
r2 = symmetric_de.rho(g2, rho_coeffs)
grid = np.linspace(0, max(g_grid)*len(rho_coeffs), r2.size//2)
sampled_de.plot_samples(r1[0], grid, n_samples)
sampled_de.plot_samples(r1[1], grid, n_samples)
plt.plot(grid, r2[0,:])
plt.plot(grid, r2[1,:])
plt.show()

#Compare gamma_inv
gi1 = sampled_de.gamma_inv(r1)
gi2 = symmetric_de.gamma_inv(r2, f_grid, g_grid)
sampled_de.plot_samples(gi1, f_grid, n_samples)
plt.plot(f_grid, gi2)
plt.show()

#Compare lambda
l1 = sampled_de.lambd(gi1, lambda_coeffs, n_samples)
l2 = symmetric_de.lambd(gi2, lambda_coeffs)
grid = np.linspace(-max(f_grid)*len(lambda_coeffs), max(f_grid)*len(lambda_coeffs), l2.size)
sampled_de.plot_samples(l1, grid, n_samples)
plt.plot(grid, l2)
plt.show()

# #Compare convolution
c1 = sampled_de.conv(l1, samples, n_samples)
c2 = symmetric_de.conv(l2, cdf)
grid = np.linspace(-max(f_grid)*6, max(f_grid)*6, c2.size)
sampled_de.plot_samples(c1, grid, n_samples)
plt.plot(grid, c2)
plt.show()