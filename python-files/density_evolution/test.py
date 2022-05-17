"""
Implementation of asymmetric density evolution from : https://arxiv.org/pdf/cs/0509014.pdf
"""
#%%
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
    g1 = 1 - interp.griddata(f_grid, f, -np.log(np.tanh(g_grid[0,:]/2)))
    g2 = interp.griddata(f_grid, f, np.log(np.tanh(g_grid[1,:]/2)))

    return np.stack((g1,g2))

def gamma_inv(g, f_grid, g_grid):
    n = f.size//2

    f1 = 1 - interp.griddata(g_grid[0,:], f[:n], g_grid[0,:], fill_value=0)
    f2 = interp.griddata(g_grid[1,:], f[n-1:], g_grid[1,:])

    return np.stack((f1, f2))

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

    return y

#initialisation
pl_0 = 0 #Cdf
pl_1 = 0 #Cdf
n_iter = 2

llr_lim = 10
n_grid_f = 2**4
n_grid_g = 2**4

#Pre calculations
f_grid, g_grid = init_grids(llr_lim, n_grid_f, n_grid_g)
#g_grid_inv = init_g_grid_inv(g_grid)

f = np.linspace(0,1, n_grid_f)
g = gamma(f, f_grid, g_grid)

plt.step(f_grid, f)
plt.show()
plt.step(g_grid[0,:], g[0,:])
plt.show()
plt.step(g_grid[1,:], g[1,:])
plt.show()
#%%
#%%
#Density evolution
for i in range(n_iter):
    #First term of (15)
    f1 = (pl_0 + pl_1)/2
    g1 = gamma(f1, f_grid, g_grid)
    rho1 = rho(g1)

    #Second term of (15)
    f2 = (pl_0 - pl_1)/2
    g2 = gamma(f2)
    rho2 = rho(g2)

    ql_0 = gamma_inv(rho1 + rho2)
    ql_1 = gamma_inv(rho1 - rho2)

    pl_0 = conv(pl_0, lamb(ql_0))
    pl_1 = conv(pl_0, lamb(ql_1))



    
