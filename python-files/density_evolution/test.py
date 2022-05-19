"""
Implementation of asymmetric density evolution from : https://arxiv.org/pdf/cs/0509014.pdf
"""
#%%
import numpy as np 
import matplotlib.pyplot as plt
import scipy.interpolate as interp
#%%
def init_grids(max_val, f_n, g_n):
    """Returns an evenly spaced grid of n points."""
    f_grid = np.append(np.linspace(-max_val, max_val, f_n-1), np.inf)

    f_step = abs(f_grid[1] - f_grid[0])
    max_val = -np.log(np.tanh(f_step))
    
    min_step = -np.log(np.tanh(max_val))
    g_n_max = int(2 ** np.floor(np.log2(max_val/min_step)))
    assert g_n < g_n_max, "Too small step size in g_grid."
    
    g_grid = np.append(np.linspace(0, max_val, g_n//2-1), np.inf)
    
    return f_grid, g_grid

def gamma(f, f_grid, g_grid):
    zero_index = f_grid.size//2 - 1

    g0 = np.array([1-f[-1]])
    g0 = np.append(g0, 1 - interp.griddata(f_grid, f, -np.log(np.tanh(g_grid/2))[1:-1]))
    g0 = np.append(g0, 1 - f[zero_index - 1])

    g1 = np.array([0])
    g1 = np.append(g1, interp.griddata(f_grid, f, np.log(np.tanh(g_grid/2))[1:-1]))
    g1 = np.append(g1, f[zero_index])

    return np.stack((g0,g1))

def gamma_inv(g, f_grid, g_grid):
    n = f_grid.size//2

    f1 = interp.griddata(g_grid, g[1,:], -np.log(np.tanh(-f_grid[:n]/2)), fill_value=0)
    f1[-1] = g[0,-1]
    f2 = interp.griddata(g_grid, g[0,:], -np.log(np.tanh(f_grid[n:]/2)),  fill_value = 0)

    return np.hstack((f1, f2))

def rho(x):

    coeffs = [0, 0, 0, 0, 0, 0.9032,  0.0968]

    x_pad = np.pad(x, [(0,0), (0, x[0,:].size)])
    x_ft = np.fft.fft(x_pad)
    
    x_pow = np.ones(x_ft.shape, dtype=np.complex)
    y_ft = np.zeros(x_ft.shape, dtype=np.complex)

    for coeff in coeffs:
        y_ft += coeff*x_pow
    
        x_pow[0,:] = x_pow[0,:] * x_ft[0,:] + x_pow[1,:] * x_ft[1,:]
        x_pow[1,:] = x_pow[1,:] * x_ft[0,:] + x_pow[0,:] * x_ft[1,:]
        
        # plt.plot(np.abs(x_pow[0,:]))
        # plt.plot(np.abs(y_ft[0,:]))
        # plt.plot(np.abs(x_ft[0,:]))
        # plt.show()
        
    
    y = np.fft.ifft(y_ft)
    y = np.abs(y)

    plt.plot(y[0,:])
    plt.plot(y[1,:])

    plt.show()

    return y

def rho(x):

    coeffs = [0, 0, 0, 0, 0, 0.9032,  0.0968]


    x_pow = np.zeros((2, x.size//2))
    x_pow[0,:] = (np.convolve(x[0,:],x[0,:]) + np.convolve(x[1,:], x[1,:]))[:x.size//2] * 20/2**11
    x_pow[1,:] = (np.convolve(x[1,:],x[0,:]) + np.convolve(x[0,:], x[1,:]))[:x.size//2] * 20/2**11
    
    plt.plot(x_pow[0,:])
    plt.plot(x_pow[1,:])
    plt.show()
    y = np.zeros(x.shape)

    for coeff in coeffs[1:]:
        y += coeff*x_pow
    
        x_pow[0,:] = (np.convolve(x_pow[0,:],x[0,:]) + np.convolve(x_pow[1,:], x[1,:]))[:x.size//2] * 20/2**12
        x_pow[1,:] = (np.convolve(x_pow[1,:],x[0,:]) + np.convolve(x_pow[0,:], x[1,:]))[:x.size//2] * 20/2**12

        plt.plot(np.abs(x_pow[0,:]))
        plt.plot(np.abs(x_pow[1,:]))
        plt.show()

        plt.plot(y[0,:])
        plt.plot(y[1,:])
        plt.title("Y")
        plt.show()
        

    plt.show()

    return y


def lambd(x):
    coeffs = [0, 0, 0.2895, 0.3158, 0, 0, 0.3947]
    
    x_ft = np.fft.fft(x)
    
    x_pow = np.ones(x.size, dtype = np.complex)
    y_ft = np.zeros(x_ft.shape, dtype = np.complex)

    for coeff in coeffs:
        y_ft += coeff*x_pow

        x_pow = x_pow * x_ft        
    
    y = np.fft.ifft(y_ft)

    return np.abs(y)

def conv(x1, x2):
    x1_ft = np.fft.fft(x1)
    x2_ft = np.fft.fft(x2)
    conv = np.fft.ifft(x1_ft*x2_ft)
    
    return np.abs(conv)

#Initialisation
llr_lim = 10
n_grid_f = 2**11
n_grid_g = 2**11
n_iter = 1

pl_0 = np.zeros(n_grid_f)
pl_0[2**10 - 20*2**3] = 1-0.91606
pl_0[2**10 + 20*2**3] = 0.91606
pl_0 = np.cumsum(pl_0)
pl_1 = np.copy(pl_0)

ql_0 = np.zeros(n_grid_f)
ql_1 = np.zeros(n_grid_f)

f_grid, g_grid = init_grids(llr_lim, n_grid_f, n_grid_g)
error = 1

# plots = [plt.subplot(3,2,i) for i in range(1,5)]
# plots.append(plt.subplot(3,2,(5,6)))

#Density evolution
for i in range(n_iter):
    
    # plots[0].step(f_grid, pl_0)
    # plots[1].step(f_grid, pl_1)
    # plots[2].step(f_grid, ql_0)
    # plots[3].step(f_grid, ql_1)
    # plots[4].scatter(i, error)

    
    #First term of (15)
    F1 = (pl_0 + pl_1)/2
    G1 = gamma(F1, f_grid, g_grid)
    rho1 = rho(G1)

#%%
    # plt.plot(f_grid, F1)
    # plt.plot(f_grid, F1)
    # plt.show()
    # plt.plot(g_grid, G1[0,:])
    # plt.plot(g_grid, G1[1,:])
    # plt.show()
    # plt.plot(g_grid, rho1[0,:])
    # plt.plot(g_grid, rho1[1,:])
    # plt.show()

    #Second term of (15)
    F2 = (pl_0 - pl_1)/2
    G2 = gamma(F2, f_grid, g_grid)
    rho2 = rho(G2)

    ql_0 = gamma_inv(rho1 + rho2, f_grid, g_grid)
    ql_1 = gamma_inv(rho1 - rho2, f_grid, g_grid)

    pl_0 = conv(pl_0, lambd(ql_0))
    pl_1 = conv(pl_1, lambd(ql_1))

    error = pl_0[f_grid == 0] + (1-pl_1[f_grid==0])


#%% Linearity check
fg, gg = init_grids(10, 2**8, 2**8)
f1 = np.linspace(0, 1, 2**8)
f2 = np.linspace(0, 1, 2**8)**2
s = 0.3
g1 = s*gamma(f1, fg, gg) + (1-s)*gamma(f2, fg, gg)
g2 = gamma(s*f1 + (1-s)*f2, fg, gg)

plt.step(fg, f1)
plt.step(fg, f2)
plt.step(gg, g1[0,:])
plt.step(gg, g2[0,:])
plt.step(gg, g1[1,:])
plt.step(gg, g2[1,:])

plt.show()
# %%
