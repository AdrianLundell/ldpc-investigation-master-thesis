"""
This file implements the full asymmetric density evolution algorithm for optimisng on the discritized A-BIAWGN channel.

Sources: 
https://arxiv.org/pdf/cs/0509014.pdf
https://arxiv.org/pdf/2001.01249.pdf

"""

import de_methods
import generate_distributions

import numpy as np 
import matplotlib.pyplot as plt

#%%
def density_evolution(p0_pdf, f_grid, g_grid, r, l, n_iter = 50, tol = 1e-6):

    # fig, axes = plt.subplots(1,2)
    p0 = symmetric_de.to_cdf(p0_pdf)
    pl = p0

    for i in range(n_iter):
        x1 = symmetric_de.gamma(pl, f_grid, g_grid)
        
        x2 = symmetric_de.rho(x1, r)
        x3 = symmetric_de.gamma_inv(x2, f_grid, g_grid)
        x4 = symmetric_de.lambd(x3, l)
        pl = symmetric_de.conv(x4, p0)

        error = pl[pl.size//2]
        # axes[0].plot(pl)
        # axes[1].scatter(i, error)

        if error < tol:
            return True

    # plt.show()
    return False

#%% Single run
n_grid = 8192
sigma, p0, bins = discretize_sigma.compute_pdf(0.01, 0.5, n_grid, 30)
f_grid, g_grid, pdf = discretize_sigma.create_pdf(p0, bins, n_grid)
cdf = symmetric_de.to_cdf(pdf)
lambda_coeffs = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])
rho_coeffs = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1])

density_evolution(pdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, 50)

#%% Many runs
import time 

t0 = time.time()
n_grid = 8192

rho_coeffs = np.array([0, 
1.1203e-04,
0.0023,
0.0056,
0.0092,
0.0126,
0.0153,
0.0172,
0.0185,
0.0193,
0.0203,
0.0220,
0.0235,
0.0253,
0.0267,
0.0251,
0.0240,
0.0195,
6.4078e-04,
0.7129])

lambda_coeffs = np.array([0, 
0.9464,
0.0325,
0.0033,
0.0022,
0.0018,
0.0015,
0.0013,
0.0012,
0.0011,
0.0010,
9.8768e-04,
9.4196e-04,
9.0299e-04,
8.6939e-04,
8.4013e-04,
8.1443e-04,
7.9168e-04,
7.7141e-04,
7.5322e-04])

for rber in np.linspace(0.001, 0.1):
    sigma, p0, bins = discretize_sigma.compute_pdf(rber, 0.5, n_grid, 30)
    f_grid, g_grid, pdf = discretize_sigma.create_pdf(p0, bins, n_grid)
    cdf = symmetric_de.to_cdf(pdf)
    
    converged = density_evolution(pdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, 50)
    
    if not converged:
        print(f"RBER : {rber}")
        break


print(f"{time.time() - t0} seconds")

# %%

def density_evolution(p0_pdf, f_grid, g_grid, n_iter = 50):

    pl = to_cdf(p0_pdf)

    for l in range(n_iter):
        plt.plot(pl)
        plt.show()
        x1 = gamma(pl, f_grid, g_grid)
        plt.plot(g_grid, x1[0,:])
        plt.show()
        plt.plot(g_grid, x1[1,:])
        plt.show()
        
        x2 = rho(x1)
        x3 = gamma_inv(x2, f_grid, g_grid)
        x4 = lambd(x3)
        pl = conv(x4, )

        plt.show()


#%% Bisection search
def bisection_search(min, max, eval, tol = 1e-5):
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