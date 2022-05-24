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
def density_to_dist(density):
    """Returns the discrete pdf of a cdf"""
    return np.cumsum(density)

def dist_to_density(dist):
    """Returns the discrete cdf of a pdf"""
    dist = np.hstack((0, dist, 1))
    density = dist[1:] - dist[:-1]
    return density[:-1]

#%% Test of probability convertion functions
N = 10
unif_cdf = np.linspace(1/N,1, N, endpoint=True)
unif_pdf = np.full(N, 1/N)

plt.plot(density_to_dist(unif_pdf))
plt.plot(unif_cdf)
plt.show()

plt.plot(dist_to_density(unif_cdf))
plt.plot(unif_pdf)
plt.show()

#%% Gamma convertion functions
def gamma(F, F_grid, G_grid):
    zero_index = F.size//2-1
    f_step = abs(F_grid[1] - F_grid[0])
    
    G0_indices = np.ceil(-np.log(np.tanh(G_grid/2)) / f_step)
    G0_indices = np.clip(G0_indices, 0, zero_index).astype(int)
    G0 = 1 - F.take(G0_indices + zero_index)

    G1_indices = np.ceil(np.log(np.tanh(G_grid/2)) / f_step)
    G1_indices = np.clip(G1_indices, -zero_index, 0).astype(int)
    G1 = F.take(G1_indices + zero_index)

    return np.stack((G0, G1))

def gamma_inv(g, f_grid, g_grid):
    n = f_grid.size//2

    f_neg = interp.griddata(g_grid, g[1,:], -np.log(np.tanh(-f_grid[:n-1]/2)))
    f_0 = g[0,-1]
    f_pos = 1-interp.griddata(g_grid, g[0,:], -np.log(np.tanh(f_grid[n:]/2)))

    return np.hstack((f_neg, f_0, f_pos))

# Gamma test
#GAMMA(distribution) is the distribution of gamma(messages)
#Simulated density and computed density should look the same
m = -10 + 20*np.random.rand(10**5)

n = 2**10
m_cdf = np.linspace(0,1, n)
F_grid = np.linspace(-10, 20, n)
G_grid = np.linspace(0, 8, n//2)
G_computed = gamma(m_cdf, F_grid, G_grid)

G0_sampled = -np.log(np.tanh(np.abs(m[m<=0])/2))
G0_hist, bins = np.histogram(G0_sampled, bins=G_grid)
G1_sampled = -np.log(np.tanh(np.abs(m[m>0])/2))
G1_hist, bins = np.histogram(G1_sampled, bins=G_grid)

plt.plot(G_grid, dist_to_density(G_computed[0,:])*10**5)
plt.plot(G_grid[:-1], G0_hist)
plt.show()
plt.plot(G_grid, dist_to_density(G_computed[1,:])*10**5, "--")
plt.plot(G_grid[:-1], G0_hist, "--")
plt.show()

#%%Gamma inverse should yield same result back
plt.plot(m_cdf)
plt.plot(gamma_inv(G, F_grid, G_grid))
plt.show()

#%% Invertion test of gamma


#%% Linearity test of gamma



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
