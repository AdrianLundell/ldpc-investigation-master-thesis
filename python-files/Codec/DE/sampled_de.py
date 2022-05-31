"""
This file contains methods for simulating density evolution by taking a large number of samples.
For comparison purposes, to slow for practical use.

Source: https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=910578
"""

import numpy as np 
import matplotlib.pyplot as plt 

def to_cdf(pdf):
    """Returns the discrete pdf of a cdf"""
    return np.cumsum(pdf)

def to_pdf(cdf):
    """Returns the discrete cdf of a pdf"""
    cdf = np.hstack((0, cdf, 1))
    pdf = cdf[1:] - cdf[:-1]
    return pdf[:-1]

def plot_samples(samples, bins, n_samples, cdf = True):
    """Compute and plot the pdf of samples"""
    bins = np.append(bins, np.inf)
    values, bins = np.histogram(samples, bins) 
    values = values/ n_samples
    
    if cdf: 
        values = to_cdf(values)
    
    plt.plot(bins[:-1], values)

def sample(cdf, f_grid, n_samples):
    """Sample n_samples from the discrete stochastic distribution determined by (cdf, f_grid)"""
    samples = []
    for i in range(n_samples):
        x = np.random.rand()
        sample = f_grid[np.argmax(cdf >= x) - 1]  
        samples.append(sample)
    return np.array(samples)

def gamma(samples):    
    g0 = -np.log(np.tanh(np.abs(samples[samples>0])/2))    
    g1 = -np.log(np.tanh(np.abs(samples[samples<=0])/2))
    return (g0, g1)

def rho(samples, rho_coeffs, n_samples):
    
    sampled_rho0 = []
    for i in range(n_samples//2):
        sample = 0
        for i, coeff in enumerate(rho_coeffs[1:]):
            x = np.random.randint(0, samples[0,:].size, i)
            x = np.take(samples[0,:], x)
            sample = coeff * np.sum(x)
        sampled_rho0.append(sample)
    rho0 = np.array(sampled_rho0)

    sampled_rho1 = []
    for i in range(n_samples//2):
        sample = 0
        for i, coeff in enumerate(rho_coeffs[1:]):
            x = np.random.randint(0, samples[1,:].size, i)
            x = np.take(samples[1,:], x)
            sample = coeff * np.sum(x)
        sampled_rho1.append(sample)
    rho1 = np.array(sampled_rho1)
    
    return np.stack((rho0, rho1))

def gamma_inv(samples):
    sampled_inv0 = -np.log((1 + np.exp(-samples[0,:]))/(1 - np.exp(-samples[0,:])))
    sampled_inv1 = np.log((1 + np.exp(-samples[1,:]))/(1 - np.exp(-samples[1,:])))
    sampled_inv = np.append(sampled_inv0, sampled_inv1)
    
    return sampled_inv

def lambd(samples, lambda_coeffs, n_samples):
    sampled_lambda = []
    for i in range(n_samples):
        sample = 0
        for i, coeff in enumerate(lambda_coeffs[1:]):
            x = np.random.randint(0, samples.size, coeff)
            x = np.take(samples, x)
            sample = coeff * np.sum(x)
        sampled_lambda.append(sample)
    
    return np.array(sampled_lambda)

def conv(samples1, samples2, n_samples):
    sampled_conv = []
    for i in range(n_samples):
        x1 = samples1[np.random.randint(0, n_samples)]
        x2 = samples2[np.random.randint(0, n_samples)]
        sample = x1 + x2
        sampled_conv.append(sample)
    sampled_conv = np.array(sampled_conv)

    return sampled_conv