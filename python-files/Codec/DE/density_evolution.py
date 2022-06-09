"""
This file implements the full asymmetric density evolution algorithm for optimisng on the discritized A-BIAWGN channel.

Sources: 
https://arxiv.org/pdf/cs/0509014.pdf
https://arxiv.org/pdf/2001.01249.pdf

"""

import time
from unicodedata import name 
import de_methods
import generate_distributions

import numpy as np 
import matplotlib.pyplot as plt
np.seterr(divide='ignore')

def symmetric_density_evolution(cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, n_iter = 50, tol = 1e-6, plot = False):
    if plot: 
        fig, axes = plt.subplots(1,2)

    #assert np.sum(rho_coeffs[1:]) == 1, "Invalid rho polynom"
    #assert np.sum(lambda_coeffs[1:]) == 1, "Invalid lambda polynom"

    pl = cdf 
    p0 = cdf
    pl_old = 0
    i = 0
    diff = np.inf
    error = np.inf
    while (diff > 1e-8) and (error > tol) and (i < n_iter):
        x1 = de_methods.gamma(pl, f_grid, g_grid)
        x2 = de_methods.rho(x1, rho_coeffs)
        x3 = de_methods.gamma_inv(x2, f_grid, g_grid)
        x4 = de_methods.lambd(x3, lambda_coeffs)
        pl = de_methods.conv(x4, p0)

        diff = sum((pl_old - pl)**2)/np.where((pl>0) & (pl <1))[0].size
        pl_old = pl

        zero_index = pl.size//2
        error = pl[zero_index]
        
        if plot:
            axes[0].plot(pl)
            axes[1].scatter(i, error)

        i += 1

    if plot:
        plt.show()

    return np.round(error, int(-np.log10(tol)))

def bisection_search(min, max, eval, tol = 1e-4):
    while max - min > tol:
        x = (min + max)/2 
        result = eval(x)

        if result == 0:
            min = x
        elif result > 0:
            max = x
        else:
            raise Exception

    return (min + max)/2

#%%
if __name__ == "__main__":
    t0 = time.time()
    n_grid = 8192
    rho_coeffs = np.array([0,1.1203e-04, 0.0023, 0.0056, 0.0092,0.0126,0.0153,0.0172,0.0185,0.0193,0.0203,0.0220,0.0235,0.0253,0.0267,0.0251,0.0240,0.0195,6.4078e-04, 0.7129])
    lambda_coeffs = np.array([0,0.9464,0.0325,0.0033,0.0022,0.0018,0.0015,0.0013,0.0012,0.0011,0.0010,9.8768e-04,9.4196e-04,9.0299e-04,8.6939e-04,8.4013e-04,8.1443e-04,7.9168e-04,7.7141e-04,7.5322e-04])

    min = 0.01
    max = 0.1

    def eval(x):
        sigma, p0, bins = generate_distributions.compute_pdf(x, 0.5, n_grid, 30)
        f_grid, g_grid, pdf = generate_distributions.create_pdf(p0, bins, n_grid)
        cdf = de_methods.to_cdf(pdf)

        return symmetric_density_evolution(cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, plot = False)

    result = bisection_search(min, max, eval)

    print(result)
    print(time.time()-t0)

#%%
