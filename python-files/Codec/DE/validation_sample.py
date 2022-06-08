"""
This file validates the implementation of the density evolution by comparing the resulting distributions to
a sampled version
"""
import de_methods_sampled as sampled_de
import de_methods as symmetric_de
import generate_distributions

import numpy as np 
import matplotlib.pyplot as plt

n_grid = 8192
sigma, p0, bins = generate_distributions.compute_pdf(0.01, 0.5, n_grid, 30)
f_grid, g_grid, pdf = generate_distributions.create_pdf(p0, bins, n_grid)
cdf = sampled_de.to_cdf(pdf)
n_samples = 10**5
samples = sampled_de.sample(cdf, f_grid, n_samples)
lambda_coeffs = np.array([0, 0, 1])
rho_coeffs = np.array([0, 0, 0, 0, 0, 1])

rho_coeffs = np.array([0,1.1203e-04, 0.0023, 0.0056, 0.0092,0.0126,0.0153,0.0172,0.0185,0.0193,0.0203,0.0220,0.0235,0.0253,0.0267,0.0251,0.0240,0.0195,6.4078e-04, 0.7129])
lambda_coeffs = np.array([0,0.9464,0.0325,0.0033,0.0022,0.0018,0.0015,0.0013,0.0012,0.0011,0.0010,9.8768e-04,9.4196e-04,9.0299e-04,8.6939e-04,8.4013e-04,8.1443e-04,7.9168e-04,7.7141e-04,7.5322e-04])


#Compare initial distributions
#sampled_de.plot_samples(samples, f_grid, n_samples)
# plt.plot(f_grid, cdf)
# plt.title("Initial distributions")
# plt.show()

#Compare gamma
g1 = sampled_de.gamma(samples)
g2 = symmetric_de.gamma(cdf, f_grid, g_grid)
# sampled_de.plot_samples(g1[0], g_grid, n_samples)
# sampled_de.plot_samples(g1[1], g_grid, n_samples)
# plt.plot(g_grid, g2[0,:])
# plt.plot(g_grid, g2[1,:])
# plt.title("Gamma")
# plt.show()

#Compare rho
r1 = sampled_de.rho(g1, rho_coeffs, n_samples)
r2 = symmetric_de.rho(g2, rho_coeffs)
# grid = np.linspace(0, max(g_grid)*len(rho_coeffs), r2.size//2)
# sampled_de.plot_samples(r1[0], grid, n_samples)
# sampled_de.plot_samples(r1[1], grid, n_samples)
# plt.plot(grid, r2[0,:])
# plt.plot(grid, r2[1,:])
# plt.title("Rho")
# plt.show()

#Compare gamma_inv
gi1 = sampled_de.gamma_inv(r1)
gi2 = symmetric_de.gamma_inv(r2, f_grid, g_grid)
# sampled_de.plot_samples(gi1, f_grid, n_samples)
# plt.plot(f_grid, gi2)
# plt.title("Gamma inverse")
# plt.show()

#Compare lambda
l1 = sampled_de.lambd(gi1, lambda_coeffs, n_samples)
import time 

t0 = time.time()
l2 = symmetric_de.lambd1(gi2, lambda_coeffs)
print(time.time() - t0)


t0 = time.time()
l2 = symmetric_de.lambd2(gi2, lambda_coeffs)
print(time.time() - t0)

#grid = np.linspace(-max(f_grid)*len(lambda_coeffs), max(f_grid)*len(lambda_coeffs), l2.size)
#sampled_de.plot_samples(l1, grid, n_samples)
plt.plot(grid, l2)
plt.title("Lambda")
plt.show()

# #Compare convolution
# c1 = sampled_de.conv(l1, samples, n_samples)
# c2 = symmetric_de.conv(l2, cdf)
# grid = np.linspace(-max(f_grid)*6, max(f_grid)*6, c2.size)
# grid = np.linspace(-max(f_grid)*4, max(f_grid)*4, c2.size)
# sampled_de.plot_samples(c1, grid, n_samples)
# plt.plot(grid, c2)
# plt.title("Convolution")
# plt.show()