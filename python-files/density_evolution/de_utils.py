"""
This file contains methods for a numeric discretized implementation of density evolution for a general distribution.

Sources: 
https://arxiv.org/pdf/cs/0509014.pdf [1]
https://arxiv.org/pdf/2001.01249.pdf [2]

"""
import matplotlib.pyplot as plt
from scipy.stats import norm
import Analysis.utils as utils
import Modem.optimal_thresholds as optimize
import numpy as np
import scipy.signal as sp
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))


global data
data = None


def to_cdf(pdf):
    """Returns the discrete pdf from a given cdf"""
    return np.cumsum(pdf)


def to_pdf(cdf):
    """Returns the discrete pdf from a given cdf"""
    cdf = np.hstack((0, cdf, 1))
    pdf = cdf[1:] - cdf[:-1]
    return pdf[:-1]


def convolution_pad(x, final_size):
    """Returns the given array padded to at least double the final_length to avoid circular convolution, rounded to the nearest power of two"""
    nearest_power = np.ceil(np.log2(final_size))
    padding = int(2**nearest_power - x.size)
    x = np.pad(x, (0, padding))

    return x


def gamma(F, F_grid, G_grid):
    """
    Given a discrete stochastic variable f defined by a cdf with probabilities F for each value on F_grid,
    calculates the cdf of the transformed variable 
        g = gamma(f)
    with gamma defined as in [1]. 
    Because G has to be an equidistant grid for numeric purposes, this creates a small quantization error.
    """
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
    """
    Given a discrete stochastic variable g defined by a 2 dimensional cdf with probabilities G for each value on (GF(2) x G_grid),
    calculates the cdf of the transformed variable 
        f = gamma^-1(g)
    with gamma^-1 defined as the inverse of gamma in [2].
    """
    zero_index = F_grid.size//2
    G_step = abs(G_grid[1] - G_grid[0])

    F_neg_indices = np.floor(-np.log(np.tanh(-F_grid[:zero_index]/2)) / G_step)
    F_neg_indices = np.clip(F_neg_indices, 0, G[1, :].size-1).astype(int)
    F_neg = G[1, :].take(F_neg_indices)

    f_0 = G[1, -1]

    F_pos_indices = np.floor(
        -np.log(np.tanh(F_grid[zero_index+1:]/2)) / G_step)
    F_pos_indices = np.clip(F_pos_indices, 0, G[0, :].size-1).astype(int)
    F_pos = 1 - G[0, :].take(F_pos_indices)

    return np.hstack((F_neg, f_0, F_pos))


def rho(x, coeffs):
    """
    """
    final_size = (x[0, :].size - 1)*len(coeffs) + 1
    dx0, dx1 = to_pdf(x[0, :]), to_pdf(x[1, :])
    x0, x1 = convolution_pad(
        x[0, :], final_size), convolution_pad(x[1, :], final_size)
    dx0, dx1 = convolution_pad(
        dx0, final_size), convolution_pad(dx1, final_size)

    x0_ft, x1_ft = np.fft.fft(x0), np.fft.fft(x1)
    dx0_ft, dx1_ft = np.fft.fft(dx0), np.fft.fft(dx1)
    y0_ft, y1_ft = np.zeros(x0.shape, complex), np.zeros(x1.shape, complex)

    for coeff in coeffs[1:]:
        x0_ft, x1_ft = x0_ft*dx0_ft + x1_ft*dx1_ft, \
            x0_ft*dx1_ft + x1_ft*dx0_ft
        y0_ft, y1_ft = y0_ft + coeff * x0_ft, y1_ft + coeff * x1_ft

    y = np.stack((np.abs(np.fft.ifft(y0_ft)[:final_size]), np.abs(
        np.fft.ifft(y1_ft)[:final_size])))
    y[0, np.argmax(y[0, :]):] = np.max(y[0, :])
    y[1, np.argmax(y[1, :]):] = np.max(y[1, :])

    return y


def lambd(x, coeffs):
    """
    """
    final_size = (x.size - 1)*len(coeffs) + 1
    dx = to_pdf(x)

    x = convolution_pad(x, final_size)
    dx = convolution_pad(dx, final_size)

    x_ft = np.fft.fft(x)
    dx_ft = np.fft.fft(dx)
    y_ft = np.zeros(x.size, complex)

    for coeff in coeffs[1:]:
        x_ft = x_ft * dx_ft
        y_ft += coeff*x_ft

    y = np.abs(np.fft.ifft(y_ft)[:final_size])
    y[np.argmax(y):] = 1
    return y


def conv(x, x0):
    """
    """
    dx = to_pdf(x)
    final_size = x.size + x0.size - 1

    x0 = convolution_pad(x0, final_size)
    dx = convolution_pad(dx, final_size)

    y = np.abs(np.fft.ifft(np.fft.fft(dx)*np.fft.fft(x0))[:final_size])
    y[np.argmax(y):] = 1
    return y


def conv(x, x0):
    """
    Naive implementation of conv using sp.convolve, kept for reference.
    """
    x0 = np.pad(x0, (x0.size, x0.size), constant_values=(0, 1))
    x = np.pad(x, (x.size, x.size), constant_values=(0, 1))
    dx = to_pdf(x)

    y = sp.convolve(x0, dx)

    current_size = y.size
    y = y[:np.argmax(y)+1]
    y = np.pad(y, (0, current_size-y.size), constant_values=y.max())

    y = y[y.size//4: -y.size//4]
    return y


def lambd(x, coeffs):
    """
    Naive implementation of lambd using sp.convolve
    """
    x = np.pad(x, (x.size//2, x.size//2), constant_values=(0, 1))
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
        x = np.pad(x, (padding1, padding2 + padding3), constant_values=(0, 1))

        y += coeff*x
        x = x[padding1:-padding2]

    y = y[y.size//4:-y.size//4]
    return y


def rho(x, coeffs):
    """

    """
    dx = np.stack((to_pdf(x[0, :]), to_pdf(x[1, :])))
    final_size = x.size//2 * len(coeffs)

    x0 = x[0, :]
    x1 = x[1, :]
    y = np.zeros((2, final_size))
    for coeff in coeffs[1:]:
        x0, x1 = sp.convolve(x0, dx[0, :]) + sp.convolve(x1, dx[1, :]),\
            sp.convolve(x1, dx[0, :]) + sp.convolve(x0, dx[1, :])
        current_size = x.size//2

        x0 = x0[:np.argmax(x0)+1]
        x1 = x1[:np.argmax(x1)+1]
        x0 = np.pad(x0, (0, final_size-x0.size), constant_values=x0.max())
        x1 = np.pad(x1, (0, final_size-x1.size), constant_values=x1.max())
        y[0, :] += coeff * x0
        y[1, :] += coeff * x1

        x0 = x0[:current_size]
        x1 = x1[:current_size]

    return y


def init_pdf(rber, n_grid=512, llr_max=30):
    """Returns density values of a DMC pdf and its grid in F and G from a look up table."""
    mu1 = -1
    mu2 = 1

    # Load database if not loaded
    global data
    if data is None:
        data = np.loadtxt(
            "C:/Users/adria/.vscode/Projects/ldpc-investigation-master-thesis/double_soft_symmetric.txt", delimiter=",", skiprows=1)

    # Interpolate values on rber
    rber_step = data[1, 1] - data[0, 1]
    rber_index = (rber - data[0, 1])/rber_step
    index0 = int(np.floor(rber_index)) - 1
    index1 = int(np.ceil(rber_index)) - 1
    values = data[index0, :] + \
        (data[index1] - data[index0]) * (rber_index - index0)

    # Set values
    sigma = values[2]
    skew = values[3]
    sigma1 = sigma * (1-skew)
    sigma2 = sigma * skew
    thresholds = values[5:8]
    llrs = values[8:]

    # Calculate probabilities
    p0 = norm.cdf(np.append(thresholds, np.inf), mu2, sigma2)
    p0 = np.ediff1d(p0, to_begin=p0[0])
    p1 = 1 - norm.cdf(np.array(thresholds), mu1, sigma1)
    p1 = np.ediff1d(p1, to_begin=p1[0])

    step = 2*llr_max / n_grid
    bins = np.round(-llrs / step) + n_grid//2
    bins = np.clip(bins, 0, n_grid-1)

    x1 = np.linspace(-llr_max, llr_max, n_grid, endpoint=False)
    y = np.zeros(x1.shape)
    np.put(y, bins.astype(int), p0)

    f_step = abs(x1[1] - x1[0])
    max_val = -np.log(np.tanh(f_step))
    x2 = np.linspace(0, max_val, n_grid)

    return x1, x2, y
