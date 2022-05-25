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
    G0_indices = np.nan_to_num(G0_indices, nan = np.inf)
    G0_indices = np.clip(G0_indices, 0, zero_index).astype(int)
    G0 = 1 - F.take(G0_indices + zero_index)

    G1_indices = np.ceil(np.log(np.tanh(G_grid/2)) / f_step)
    G1_indices = np.nan_to_num(G0_indices, nan = -np.inf)
    G1_indices = np.clip(G1_indices, -zero_index, 0).astype(int)
    G1 = F.take(G1_indices + zero_index)

    return np.stack((G0, G1))

def gamma_inv(g, f_grid, g_grid):
    n = f_grid.size//2

    f_neg = interp.griddata(g_grid, g[1,:], -np.log(np.tanh(-f_grid[:n-1]/2)))
    f_0 = g[0,-1]
    f_pos = 1-interp.griddata(g_grid, g[0,:], -np.log(np.tanh(f_grid[n:]/2)))

    return np.hstack((f_neg, f_0, f_pos))


#%% Test density evolution
def plot_samples(samples, bins, ax, cdf = False):
    bins = np.append(bins, np.inf)
    values, bins = np.histogram(samples, bins) 
    values = values / len(samples)
    if cdf: 
        values = density_to_dist(values)
    ax.plot(bins[:-1], values)

def test(cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, n_samples = 10**4):

    #compare sampled and theoretical cdf
    samples = []
    for i in range(n_samples):
        x = np.random.rand()
        sample = f_grid[np.argmax(cdf >= x)]  
        samples.append(sample)
    samples = np.array(samples)

    fig, ax = plt.subplots(1,1)
    ax.set_title("Compare distributions")
    plot_samples(samples, f_grid, ax, cdf = True)
    plt.plot(f_grid, cdf)
    plt.show()

    #Compare sampled and theoretical gamma(cdf)
    sampled_g0 = -np.log(np.tanh(np.abs(samples[samples<=0])/2))
    sampled_g1 = -np.log(np.tanh(np.abs(samples[samples>0])/2))    
    g = gamma(cdf, f_grid, g_grid)
    
    fig, axes = plt.subplots(1,2, figsize = (15,10))
    fig.suptitle("Compare gamma(distribution)")
    plot_samples(sampled_g0, g_grid, axes[0])
    plot_samples(sampled_g1, g_grid, axes[1])
    axes[0].plot(g_grid, dist_to_density(g[0,:]))
    axes[1].plot(g_grid, dist_to_density(g[1,:]))
    plt.show()

    #Compare sampled and theoretical rho(gamma(cdf)). 
    sampled_rho0 = []
    for i in range(n_samples//2):
        sample = 0
        for i, coeff in enumerate(rho_coeffs[1:]):
            x = np.random.randint(0, sampled_g0.size, i)
            x = np.take(sampled_g0, x)
            sample = coeff * np.sum(x)
        sampled_rho0.append(sample)
    sampled_rho0 = np.array(sampled_rho0)

    sampled_rho1 = []
    for i in range(n_samples//2):
        sample = 0
        for i, coeff in enumerate(rho_coeffs[1:]):
            x = np.random.randint(0, sampled_g1.size, i)
            x = np.take(sampled_g1, x)
            sample = coeff * np.sum(x)
        sampled_rho1.append(sample)
    sampled_rho1 = np.array(sampled_rho1)
    
    fig, axes = plt.subplots(1,2, figsize = (15,10))
    fig.suptitle("Compare rho(gamma(distribution))")
    plot_samples(sampled_rho0, g_grid, axes[0])
    plot_samples(sampled_rho1, g_grid, axes[1])
    plt.show()

    #Compare sampled and theoretical inv gamma
    sampled_inv0 = -np.log((1 + np.exp(-sampled_rho0))/(1 - np.exp(-sampled_rho0)))
    sampled_inv1 = np.log((1 + np.exp(-sampled_rho1))/(1 - np.exp(-sampled_rho1)))
    sampled_inv = np.append(sampled_inv0, sampled_inv1)
    # sampled_pdf, bins = np.histogram(sampled_inv, 30)
    
    # plt.plot(bins[:-1], density_to_dist(sampled_pdf))
    # plt.title("Compare inv_gamma(rho(gamma(distribution)))")
    # plt.show()
    
    #Compare lambda
    sampled_lambda = []
    for i in range(n_samples):
        sample = 0
        for i, coeff in enumerate(lambda_coeffs[1:]):
            x = np.random.randint(0, sampled_inv.size, coeff)
            x = np.take(sampled_inv, x)
            sample = coeff * np.sum(x)
        sampled_lambda.append(sample)
    sampled_lambda = np.array(sampled_lambda)

    #Compare convolution
    sampled_conv = []
    for i in range(n_samples):
        x1 = sampled_inv[np.random.randint(0, n_samples)]
        x2 = samples[np.random.randint(0, n_samples)]
        sample = x1 + x2
        sampled_conv.append(sample)
    sampled_conv = np.array(sampled_conv)

    error = sum(sampled_conv == 0)/n_samples
    print("Error ", error)

n = 2**8
rho_coeffs = np.array([0, 0, 0, 0, 0, 1])
lambda_coeffs = np.array([0, 0, 1])
m_cdf = np.linspace(0.4, 1, 2)
#f_grid = np.linspace(-10, 10, n)
#g_grid = np.linspace(0, 8, n//2)

f_grid = np.array([0, np.inf])
g_grid = np.array([0, np.inf])


test(m_cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs)

#%% Invertion test of gamma


#%% Linearity test of gamma



def rho(x, coeffs):
    dx = np.stack((dist_to_density(x[0,:]), dist_to_density(x[1,:])))
    
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

def lambd(x, coeffs):
    dx = dist_to_density(x)

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
