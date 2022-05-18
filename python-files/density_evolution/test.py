"""
Implementation of asymmetric density evolution from : https://arxiv.org/pdf/cs/0509014.pdf
"""
#%%
from cProfile import label
import numpy as np 
import matplotlib.pyplot as plt
import scipy.interpolate as interp
#%%

def init_grids(max_val, f_n, g_n):
    """Returns an evenly spaced grid of n points """
    f_grid = np.append(np.linspace(-max_val, max_val, f_n-1), np.inf)

    f_step = abs(f_grid[1] - f_grid[0])
    max_val = -np.log(np.tanh(f_step))
    g_grid = np.append(np.linspace(0, max_val, g_n//2-1), np.inf)
    g_grid = np.stack((g_grid, g_grid)) 
    
    return f_grid, g_grid

def gamma(f, f_grid, g_grid):
    g1 = 1 - interp.griddata(f_grid, f, -np.log(np.tanh(g_grid[0,:]/2)), fill_value=0)
    g2 = interp.griddata(f_grid, f, np.log(np.tanh(g_grid[1,:]/2)), fill_value=0)

    return np.stack((g1,g2))

def gamma_inv(g, f_grid, g_grid):
    n = f_grid.size//2

    f1 = interp.griddata(g_grid[0,:], g[0,:], -np.log(np.tanh(-f_grid[:n]/2)))
    f1[-1] = g[0,-1]
    f2 = interp.griddata(g_grid[1,:], g[0,:], -np.log(np.tanh(f_grid[n:]/2)))

    return np.hstack((f1, f2))

def rho(x):
    coeffs = [0, 0.5, 0.5]
    x_ft = np.fft.fft(x)
    
    x_pow = np.ones((2, x.size//2))
    y_ft = np.zeros(x_ft.shape)

    for coeff in coeffs:
        x_pow[0,:] = x_pow[0,:] * x[0,:] + x_pow[1,:] * x[1,:]
        x_pow[1,:] = x_pow[1,:] * x[0,:] + x_pow[0,:] * x[1,:]
        
        y_ft += coeff*x_pow
    
    y = np.fft.ifft(y_ft)
    y = np.abs(y)

    return y

def lambd(x):
    coeffs = [0, 0.5, 0.5]
    x_ft = np.fft.fft(x)
    
    x_pow = np.ones(x.size)
    y_ft = np.zeros(x_ft.shape)

    for coeff in coeffs:
        x_pow = x_pow * x        
        y_ft += coeff*x_pow
    
    y = np.fft.ifft(y_ft)

    return np.abs(y)

def conv(x1, x2):
    x1_ft = np.fft.fft(x1)
    x2_ft = np.fft.fft(x2)
    conv = np.fft.ifft(x1_ft*x2_ft)
    
    return np.abs(conv)

#initialisation
llr_lim = 10
n_grid_f = 2**4
n_grid_g = 2**4
n_iter = 4

pl_0 = np.linspace(0,1, n_grid_f)
pl_1 = np.linspace(0,1, n_grid_f)
ql_0 = np.zeros(n_grid_f)
ql_1 = np.zeros(n_grid_f)

f_grid, g_grid = init_grids(llr_lim, n_grid_f, n_grid_g)
error = 0

plots = [plt.subplot(3,2,i) for i in range(1,5)]
plots.append(plt.subplot(3,2,(5,6)))

#Density evolution
for i in range(n_iter):
    
    plots[0].step(f_grid, pl_0)
    plots[1].step(f_grid, pl_1)
    plots[2].step(f_grid, ql_0)
    plots[3].step(f_grid, ql_1)
    plots[4].scatter(i, error)
    
    #First term of (15)
    f1 = (pl_0 + pl_1)/2
    g1 = gamma(f1, f_grid, g_grid)
    rho1 = rho(g1)

    #Second term of (15)
    f2 = (pl_0 - pl_1)/2
    g2 = gamma(f2, f_grid, g_grid)
    rho2 = rho(g2)

    ql_0 = gamma_inv(rho1 + rho2, f_grid, g_grid)
    ql_1 = gamma_inv(rho1 - rho2, f_grid, g_grid)

    pl_0 = conv(pl_0, lambd(ql_0))
    pl_1 = conv(pl_0, lambd(ql_1))

    error = sum(pl_0[f_grid < 0]) + sum(pl_1[f_grid > 0])

    

# %%
