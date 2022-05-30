#%%
import sys 
import os 
sys.path.insert(1, os.path.join(sys.path[0], '../..'))

import Modem.optimal_thresholds as optimize
import Analysis.utils as utils

import numpy as np 
import matplotlib.pyplot as plt 

#%%

def init_pdf(n_thresholds, n_grid):
    """Returns density values of a DMC pdf and its grid in F and G from a look up table."""
    pass


#%%

def compute_pdf(sigma, ratio, n_thresholds, n_grid, llr_max):
    """Returns values and bin numbers of the pdf defining the A-BIAWGN (sigma, ratio) 
       discretized by n_thresholds and quantized to the grid of n_grid points over (-llr_max, llr_max)"""

    thresholds = optimize.optimize_threhsolds(sigma)
    llrs = optimize.LLR(thresholds)
    real_sigma = utils.rber_to_sigma(sigma, ratio)

    step = 2*llr_max / n_grid
    bins = np.round(thresholds / step)

    return sigma, llrs, bins