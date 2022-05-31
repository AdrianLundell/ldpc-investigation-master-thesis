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

sigma, p0, bins = discretize_sigma.compute_pdf(0.01, 0.5, 1024, 10)
f_grid, g_grid, pdf = discretize_sigma.create_pdf(p0, bins, 1024)
cdf = sampled_de.to_cdf(pdf)
n_samples = 10**5
samples = sampled_de.sample(cdf, f_grid, n_samples)

#Compare initial distributions
# sampled_de.plot_samples(samples, f_grid, n_samples)
# plt.plot(f_grid, cdf)
# plt.show()

#Compare gamma
g1 = sampled_de.gamma(samples)
sampled_de.plot_samples(g1[0], g_grid, n_samples)
sampled_de.plot_samples(g1[1], g_grid, n_samples)
g2 = symmetric_de.gamma(cdf, f_grid, g_grid)
plt.plot(g_grid, g2[0,:])
plt.plot(g_grid, g2[1,:])
plt.show()
