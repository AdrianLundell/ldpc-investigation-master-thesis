"""
Implementation of asymmetric density evolution from : https://arxiv.org/pdf/cs/0509014.pdf
"""
#%%
import numpy as np 
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.signal as sp 
import scipy.stats as stats

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
def gamma(F, F_grid, G_grid):
    zero_index = F_grid.size//2-1
    F_step = abs(F_grid[1] - F_grid[0])
    
    G0_indices = np.floor(-np.log(np.tanh(G_grid/2)) / F_step)
    G0_indices = np.clip(G0_indices, 0, zero_index).astype(int)
    G0 = 1 - F.take(G0_indices + zero_index)

    G1_indices = np.floor(np.log(np.tanh(G_grid/2)) / F_step)
    G1_indices = np.clip(G1_indices, -zero_index, -1).astype(int)
    G1 = F.take(G1_indices + zero_index)

    return np.stack((G0, G1))

def gamma_inv(G, F_grid, G_grid):
    zero_index = F_grid.size//2-1
    G_step = abs(G_grid[1] - G_grid[0])

    F_neg_indices = np.floor(-np.log(np.tanh(-F_grid[:zero_index]/2)) / G_step)
    F_neg_indices = np.clip(F_neg_indices, 0, G[1,:].size-1).astype(int)
    F_neg = G[1,:].take(F_neg_indices)

    f_0 = G[0,-1]
    
    F_pos_indices = np.floor(-np.log(np.tanh(F_grid[zero_index+1:]/2)) / G_step)
    F_pos_indices = np.clip(F_pos_indices, 0, G[0,:].size-1).astype(int)
    F_pos = G[0,:].take(F_pos_indices)

    return np.hstack((F_neg, f_0, F_pos))

#%% Convolution functions
def rho(x):
    dx = np.stack((to_pdf(x[0,:]), to_pdf(x[1,:])))
    coeffs = np.array([0, 0, 0, 0, 0, 1])
    
    n = x.size//2
    next_power = np.ceil(np.log2(n*len(coeffs)))
    cdf_pad_size = int(2**next_power - n)
    zero_pad_size = int(2**(next_power))*3

    last_value = x[:,-1]
    x = np.pad(x, [(0,0),(0, zero_pad_size + cdf_pad_size)], constant_values = 0)
    x[0,n:zero_pad_size] = last_value[0]
    x[1,n:zero_pad_size] = last_value[1]
    dx = np.pad(dx, [(0,0),(0, zero_pad_size + cdf_pad_size)], constant_values = 0)

    x_ft = np.fft.fft(x)
    y_ft = np.zeros(x_ft.shape, dtype=complex)
    dx_ft = np.fft.fft(dx)

    # plt.plot(x[1,:])
    # plt.show()
    # plt.plot(dx[1,:])
    # plt.show()

    # plt.plot(x_ft[1,:])
    # plt.show()
    # plt.plot(dx_ft[1,:])
    # plt.show()
    
    for i, coeff in enumerate(coeffs[1:]):
        x_ft[0,:] = (x_ft[0,:]*dx_ft[0,:] + x_ft[1,:]*dx_ft[1,:])
        x_ft[1,:] = (x_ft[1,:]*dx_ft[0,:] + x_ft[0,:]*dx_ft[1,:])
        plt.plot(np.abs(x_ft[1,:]))
        plt.show()
        plt.plot(np.abs(np.fft.ifft(x_ft[1,:])))
        plt.title(i)
        plt.show()
        
        y_ft += coeff*x_ft

    y = np.abs(np.fft.ifft(y_ft))[:, :zero_pad_size/4]
    return y

def lambd(x):
    dx = to_pdf(x)
    coeffs = np.array([0, 0, 1])

    n = x.size
    next_power = np.ceil(np.log2(n*len(coeffs)))
    cdf_pad_size = int(2**next_power - n)
    zero_pad_size = int(2**next_power)

    x = np.pad(x, [(0, cdf_pad_size + zero_pad_size)], constant_values = 0)
    x[n:zero_pad_size] = x.max()
    dx = np.pad(dx, [(0, zero_pad_size + cdf_pad_size)], constant_values = 0)

    x_ft = np.fft.fft(x)
    y_ft = np.zeros(x_ft.shape, dtype=complex)
    dx_ft = np.fft.fft(dx)

    # plt.plot(x)
    # plt.show()
    # plt.plot(dx)
    # plt.show()

    # plt.plot(x_ft)
    # plt.show()
    # plt.plot(dx_ft)
    # plt.show()

    for coeff in coeffs[1:]:
        x_ft = x_ft * dx_ft
        y_ft += coeff*x_ft

    y = np.abs(np.fft.ifft(y_ft))
    return y[:n*len(coeffs)]  

#%% Run simulation
#Initialisation
def init_grids(max_val, f_n, g_n):
    """Returns an evenly spaced grid of n points."""
    f_grid = np.linspace(-max_val, max_val, f_n)

    f_step = abs(f_grid[1] - f_grid[0])
    max_val = -np.log(np.tanh(f_step))

    g_grid = np.linspace(0, max_val, g_n//2)

    return f_grid, g_grid

llr_lim = 10
n_grid_f = 512
n_grid_g = 512
n_iter = 2

#Hard read of all zero codeword
p0 = np.zeros(n_grid_f)
p0[256 - 64] = 0.3
p0[256 + 64] = 0.7


p0 = to_cdf(p0)
pl = np.copy(p0)
f_grid, g_grid = init_grids(llr_lim, n_grid_f, n_grid_g)

for l in range(n_iter):
    x1 = gamma(pl, f_grid, g_grid)
    x2 = rho(x1)
    x3 = gamma_inv(x2, f_grid, g_grid)
    x4 = lambd(x3)
    pl = sp.convolve(p0, x4)

#    plt.plot(f_grid, pl[])
# %%
