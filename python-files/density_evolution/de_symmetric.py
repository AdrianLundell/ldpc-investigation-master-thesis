"""
Implementation of asymmetric density evolution from : https://arxiv.org/pdf/cs/0509014.pdf
"""
#%%
import numpy as np 
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.signal as sp 

#%%
def density_to_dist(density):
    return np.cumsum(density)

def dist_to_density(dist):
    dist = np.hstack((0, dist, 1))
    density = dist[1:] - dist[:-1]
    return density[:-1]

def init_grids(max_val, f_n, g_n):
    """Returns an evenly spaced grid of n points."""
    f_grid = np.append(np.linspace(-max_val, max_val, f_n-1), np.inf)

    f_step = abs(f_grid[1] - f_grid[0])
    max_val = -np.log(np.tanh(f_step))
    
    min_step = -np.log(np.tanh(max_val))
    g_n_max = int(2 ** np.floor(np.log2(max_val/min_step)))
    #assert g_n < g_n_max, "Too small step size in g_grid."
    
    g_grid = np.append(np.linspace(0, max_val, g_n//2-1), np.inf)
    
    return f_grid, g_grid

def gamma(f, f_grid, g_grid):
    zero_index = f_grid.size//2 - 1

    g0 = np.array([1-f[-1]])
    g0 = np.append(g0, 1 - interp.griddata(f_grid, f, -np.log(np.tanh(g_grid[1:-1]/2))))
    g0 = np.append(g0, 1 - f[zero_index - 1])

    g1 = np.array([0])
    g1 = np.append(g1, interp.griddata(f_grid, f, np.log(np.tanh(g_grid[1:-1]/2))))
    g1 = np.append(g1, f[zero_index])

    return np.stack((g0,g1))

def gamma_inv(g, f_grid, g_grid):
    n = f_grid.size//2

    f_neg = interp.griddata(g_grid, g[1,:], -np.log(np.tanh(-f_grid[:n-1]/2)))
    f_0 = g[0,-1]
    f_pos = 1-interp.griddata(g_grid, g[0,:], -np.log(np.tanh(f_grid[n:]/2)))

    return np.hstack((f_neg, f_0, f_pos))


def rho(x):
    dx = np.stack((dist_to_density(x[0,:]), dist_to_density(x[1,:])))
    coeffs = np.array([0, 0, 0, 0, 0, 1])
    
    n = x.size//2
    zero_pad_size = n*2**len(coeffs)
    cdf_pad_size = zero_pad_size - n 

    x = np.pad(x, [(0,0),(0, cdf_pad_size)], constant_values = x[:,-1])
    x = np.pad(x, [(0,0),(0, zero_pad_size)], constant_values = 0)
    dx = np.pad(dx, [(0,0),(0, zero_pad_size + cdf_pad_size)], constant_values = 0)

    x_ft = np.fft.fft(x)
    y_ft = np.zeros(x_ft.shape, dtype=np.complex)
    dx_ft = np.fft.fft(dx)
    
    for coeff in coeffs[1:]:
        x_ft[0,:] = (x_ft[0,:]*dx_ft[0,:] + x_ft[1,:]*dx_ft[1,:])/2
        x_ft[1,:] = (x_ft[1,:]*dx_ft[0,:] + x_ft[0,:]*dx_ft[1,:])/2
        
        y_ft += coeff*x_ft

    y = np.abs(np.fft.ifft(y_ft))
    return y[:cdf_pad_size]

def lambd(x):
    dx = dist_to_density(x)
    coeffs = np.array([0, 0, 1])

    n = x.size
    zero_pad_size = n*2**len(coeffs)
    cdf_pad_size = zero_pad_size - n 

    x = np.pad(x, [(0, cdf_pad_size)], constant_values = x[-1])
    x = np.pad(x, [(0, zero_pad_size)], constant_values = 0)
    dx = np.pad(dx, [(0, zero_pad_size + cdf_pad_size)], constant_values = 0)

    x_ft = np.fft.fft(x)
    y_ft = np.zeros(x_ft.shape, dtype=np.complex)
    dx_ft = np.fft.fft(dx)

    for coeff in coeffs[1:]:
        x_ft = x_ft * dx_ft
        y_ft += coeff*x_ft

    y = np.abs(np.fft.ifft(y_ft))
    return y[:cdf_pad_size]  

#Initialisation
llr_lim = 10
n_grid_f = 2
n_grid_g = 2
n_iter = 10
epsilon = 0.9

p0 = np.zeros(n_grid_f)
p0[1] = epsilon

p0 = density_to_dist(p0)
pl = np.copy(p0)

f_grid, g_grid = init_grids(llr_lim, n_grid_f, n_grid_g)

for l in range(n_iter):
    x1 = gamma(pl, f_grid, g_grid)
    x2 = rho(x1)
    x3 = gamma_inv(x2, f_grid, g_grid)
    x4 = lambd(x3)
    pl = sp.convolve(p0, x4)

    plt.plot(x2)
# %%
