"""
This file contains methods for a numeric discretized implementation of density evolution for a general distribution.

Sources: 
https://arxiv.org/pdf/cs/0509014.pdf
https://arxiv.org/pdf/2001.01249.pdf

"""
import numpy as np
import scipy.signal as sp

def to_cdf(pdf):
    """Returns the discrete pdf of a cdf"""
    return np.cumsum(pdf)

def to_pdf(cdf):
    """Returns the discrete cdf of a pdf"""
    cdf = np.hstack((0, cdf, 1))
    pdf = cdf[1:] - cdf[:-1]
    return pdf[:-1]

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

def lambd_regular(x, coeffs):
    x = np.pad(x, (x.size//2, x.size//2), constant_values = (0,1))
    dx = to_pdf(x)
    final_size = x.size * len(coeffs)

    y = np.zeros(final_size)
    for coeff in coeffs[1:]:
        x = sp.convolve(x, dx)
        current_size = x.size

        x = x[:np.argmax(x)+1]
        padding1 = int(np.ceil((-current_size+final_size)/2))
        padding2 = int(np.floor((-current_size+final_size)/2))
        padding3 = int((current_size-x.size))
        x = np.pad(x, (padding1, padding2 + padding3), constant_values = (0,1))

        y += coeff*x
        x = x[padding1:-padding2]
    
    y = y[y.size//4:-y.size//4]
    return y

def lambd(x, coeffs):
    dx = to_pdf(x)
    final_size = x.size
    for i in range(len(coeffs)-1):
        final_size = (final_size + x.size) - 1
    
    x = np.pad(x, (0, 2**16 - x.size))
    dx = np.pad(dx, (0, 2**16 - dx.size))

    x_ft = np.fft.fft(x)
    dx_ft = np.fft.fft(dx)
    y_ft = np.zeros(2**16, complex)

    for coeff in coeffs[1:]:
        x_ft = x_ft * dx_ft
        y_ft += coeff*x_ft
    
    y = np.abs(np.fft.ifft(y_ft)[:final_size])
    y[np.argmax(y):] = 1

    return y

def conv(x, x0):
    dx = to_pdf(x)
    final_size = x.size + x0.size - 1

    x0 = np.pad(x0, (0, 2**16-x0.size))
    dx = np.pad(dx, (0, 2**16-dx.size))

    y = np.abs(np.fft.ifft(np.fft.fft(dx)*np.fft.fft(x0)))    
    y = y[:final_size]
    y[np.argmax(y):] = 1

    return y