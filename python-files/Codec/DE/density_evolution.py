#%%
import symmetric_de
import discretize_sigma

import numpy as np 
import matplotlib.pyplot as plt

#%%
def density_evolution(p0_pdf, f_grid, g_grid, r, l, n_iter = 50, tol = 1e-6):

    #fig, axes = plt.subplots(1,2)
    p0 = symmetric_de.to_cdf(p0_pdf)
    pl = p0

    for i in range(n_iter):
        x1 = symmetric_de.gamma(pl, f_grid, g_grid)
        
        x2 = symmetric_de.rho(x1, r)
        x3 = symmetric_de.gamma_inv(x2, f_grid, g_grid)
        x4 = symmetric_de.lambd(x3, l)
        pl = symmetric_de.conv(x4, p0)

        error = pl[pl.size//2]
        #axes[0].plot(pl)
        #axes[1].scatter(i, error)

        if error < tol:
            return True

   # plt.show()
    return False

#%% Single run
n_grid = 8192
sigma, p0, bins = discretize_sigma.compute_pdf(0.1, 0.5, n_grid, 30)
f_grid, g_grid, pdf = discretize_sigma.create_pdf(p0, bins, n_grid)
cdf = symmetric_de.to_cdf(pdf)
lambda_coeffs = np.array([0, 0, 1])
rho_coeffs = np.array([0, 0, 0, 0, 0, 1])

density_evolution(pdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, 50)

#%% Many runs
n_grid = 8192
rho_coeffs = np.array([0, 0, 0, 0, 0, 1])
lambda_coeffs = np.array([0, 0, 1])

for rber in np.linspace(0.001, 0.1):
    sigma, p0, bins = discretize_sigma.compute_pdf(rber, 0.5, n_grid, 30)
    f_grid, g_grid, pdf = discretize_sigma.create_pdf(p0, bins, n_grid)
    cdf = symmetric_de.to_cdf(pdf)
    
    converged = density_evolution(pdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, 50)
    
    if not converged:
        print(f"RBER : {rber}")
        break

# %%
