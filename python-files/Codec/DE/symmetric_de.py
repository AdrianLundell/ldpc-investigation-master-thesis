"""
This file contains methods for a numeric discretized implementation of a general symmetric density evolution.

Source : https://arxiv.org/pdf/cs/0509014.pdf
"""
#%%
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sp
import discretize_sigma as d

#%% Probability convertion functions
def to_cdf(pdf):
    """Returns the discrete pdf of a cdf"""
    return np.cumsum(pdf)

def to_pdf(cdf):
    """Returns the discrete cdf of a pdf"""
    cdf = np.hstack((0, cdf, 1))
    pdf = cdf[1:] - cdf[:-1]
    return pdf[:-1]

#%% Gamma convertion functions

#%%
def gamma(F, F_grid, G_grid):
    zero_index = F.size//2
    F_step = abs(F_grid[1] - F_grid[0])

    G0_indices = np.floor(-np.log(np.tanh(G_grid/2)) / F_step)
    G0_indices = np.clip(G0_indices, 0, zero_index - 1).astype(int)
    G0 = 1 - F.take(G0_indices + zero_index)

    G1_indices = np.floor(np.log(np.tanh(G_grid/2)) / F_step)
    G1_indices = np.clip(G1_indices, -zero_index, -1).astype(int)
    G1 = F.take(G1_indices + zero_index)

    return np.stack((G0, G1))

def gamma_inv(G, F_grid, G_grid):
    zero_index = F_grid.size//2
    G_step = abs(G_grid[1] - G_grid[0])

    F_neg_indices = np.floor(-np.log(np.tanh(-F_grid[:zero_index]/2)) / G_step)
    F_neg_indices = np.clip(F_neg_indices, 0, G[1,:].size-1).astype(int)
    F_neg = G[1,:].take(F_neg_indices)

    f_0 = G[1,-1]

    F_pos_indices = np.floor(-np.log(np.tanh(F_grid[zero_index+1:]/2)) / G_step)
    F_pos_indices = np.clip(F_pos_indices, 0, G[0,:].size-1).astype(int)
    F_pos = 1 - G[0,:].take(F_pos_indices)

    return np.hstack((F_neg, f_0, F_pos))

#%% Convolution functions
def rho(x, coeffs):
    dx = np.stack((to_pdf(x[0,:]), to_pdf(x[1,:])))
    final_size = x.size//2 * len(coeffs)

    x0 = x[0,:]
    x1 = x[1,:]
    y = np.zeros((2,final_size))
    for coeff in coeffs[1:]:
        x0, x1 = sp.convolve(x0, dx[0,:]) + sp.convolve(x1, dx[1,:]),\
                 sp.convolve(x1, dx[0,:]) + sp.convolve(x0, dx[1,:])
        current_size = x.size//2

        x0 = x0[:np.argmax(x0)+1]
        x1 = x1[:np.argmax(x1)+1]
        x0 = np.pad(x0, (0,final_size-x0.size), constant_values = x0.max())
        x1 = np.pad(x1, (0,final_size-x1.size), constant_values = x1.max())
        y[0,:] += coeff * x0
        y[1,:] += coeff * x1

        x0 = x0[:current_size]
        x1 = x1[:current_size]

    return y

def lambd(x, coeffs):
    x = np.pad(x, (0,x.size), constant_values = 1)
    dx = to_pdf(x)
    final_size = x.size * len(coeffs)

    # plt.plot(dx)
    # plt.show()
    # plt.plot(x)
    # plt.show()

    y = np.zeros(final_size)
    for coeff in coeffs[1:]:
        x = sp.convolve(x, dx)
        current_size = x.size

        plt.plot(x)
        plt.show()

        x = x[:np.argmax(x)+1]
        x = np.pad(x, (0,final_size-x.size), constant_values = x.max())

        y += coeff*x
        zero = x.size//2
        x = x[:current_size + 1]

        plt.plot(x)
        plt.show()

    return y

#%% Run simulation
#Initialisation
def init_grids(max_val, f_n, g_n):
    """Returns an evenly spaced grid of n points"""
    f_grid = np.linspace(-max_val, max_val, f_n)

    f_step = abs(f_grid[1] - f_grid[0])
    max_val = -np.log(np.tanh(f_step))

    g_grid = np.linspace(0, max_val, g_n//2)

    return f_grid, g_grid


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

        pl = sp.convolve(p0_pdf, x4)
        current_size = pl.size
        pl = pl[:np.argmax(pl)+1]
        pl = np.pad(pl, (0, current_size-pl.size), constant_values = pl.max())

        #plt.plot(f_grid, pl[512-128:512+128])
        #axes[1].scatter(l, pl[(2**8*4)//2-1])

        # plt.show()
        # plt.plot(x2[0,:])
        # plt.show()
        # plt.plot(x2[1,:])
        # plt.show()
        # plt.plot(x3)
        # plt.show()

        plt.show()



