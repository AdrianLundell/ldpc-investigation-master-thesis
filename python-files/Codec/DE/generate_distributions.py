"""
This file contains methods for generating the pdf of LLRS for an optimally discretized A-BIAWGN channel.

Note that the distribution is first discretized and then quantized to an equidistant grid, meaning the
final distribution is approximate. 
"""
#Path hack
import sys 
import os 
sys.path.insert(1, os.path.join(sys.path[0], '../..'))

import Modem.optimal_thresholds as optimize
import Analysis.utils as utils

import numpy as np 
from scipy.stats import norm

def init_pdf(n_thresholds, n_grid):
    """Returns density values of a DMC pdf and its grid in F and G from a look up table."""
    pass

def create_pdf(p0, bins, n_grid = 512, llr_max = 10):
    """Creates a pdf and grid from values"""
    x1 = np.linspace(-llr_max, llr_max, n_grid, endpoint = False)
    y = np.zeros(x1.shape)
    np.put(y, bins.astype(int), p0)

    f_step = abs(x1[1] - x1[0])
    max_val = -np.log(np.tanh(f_step))
    x2 = np.linspace(0, max_val, n_grid)

    return x1, x2, y    

def compute_pdf(rber, skew = 0.5, n_grid = 512, llr_max = 30, mu1 = -1, mu2 = 1):
    """Returns values and bin numbers of the pdf defining the A-BIAWGN (sigma, ratio) 
       discretized by three thresholds and quantized to the grid of n_grid points over (-llr_max, llr_max)"""
    sigma = utils.rber_to_sigma(rber, skew)
    sigma1 = (1-skew) * sigma
    sigma2 = skew * sigma
    
    thresholds, max_mi = optimize.optimize_thresholds(sigma1, sigma2, mu1, mu2, symmetric = True)
    #thresholds = [optimize.mid_point(sigma1, sigma2, mu1, mu2)]

    p0 = norm.cdf(np.append(thresholds, np.inf), mu2, sigma2)
    p0 = np.ediff1d(p0, to_begin = p0[0])
    p1 = 1 - norm.cdf(np.array(thresholds), mu1, sigma1)
    p1 = np.ediff1d(p1, to_begin = p1[0])

    llrs = optimize.llrs(thresholds, sigma1, sigma2)
    step = 2*llr_max / n_grid
    bins = np.round(-llrs / step) + n_grid//2
    bins = np.clip(bins, 0, n_grid-1)

    return sigma, p0, bins

#Create database over thresholds
if __name__ == "__main__":

    name = "double_soft_symmetric.txt"
    skew = 0.5
    rber_list = np.linspace(0.001, 0.5, 1000, endpoint=False)
    result = np.zeros((rber_list.size, 10))

    for i, rber in enumerate(rber_list):
        sigma = utils.rber_to_sigma(rber, skew)
        sigma1 = (1-skew) * sigma
        sigma2 = skew * sigma
        thresholds, capacity = optimize.optimize_thresholds(sigma1, sigma2, symmetric = True)
        llrs = optimize.llrs(thresholds, sigma1, sigma2)

        row = np.hstack((sigma, skew, capacity, thresholds, llrs))

        result[i,:] = row
    
    np.savetxt()