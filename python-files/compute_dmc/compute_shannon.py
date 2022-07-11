"""
This scripts finds the maximum decodeable noise level given a code rate for three different channels:
1. BSC
2. Discretized AWGN [only for small code rates]
3. AWGN
"""
#%% Imports
import dmc_utils
import numpy as np 
import scipy.stats as stats 

# %% Binary symmetric channel
code_rate = 1 - 9/73

min = 0.001
max = 0.5
tol = 1e-7

while max - min > tol:
    rber = (min + max)/2

    capacity = 1 - (rber * np.log2(1/rber) + (1-rber)*np.log2(1/(1-rber)))
    result = capacity - code_rate

    if result > 0:
        min = rber
    elif result < 0:
        max = rber
    else:
        raise Exception

rber = (min + max)/2
sigma_tot = dmc_utils.rber_to_sigma(rber, 0.5)
print(f"The shannon limit for the BSC with code rate ({code_rate:.4}) is {rber:.4}(RBER)/{sigma_tot:.4}(SIGMA)")

#%% Discretized AWGN
###########################################
####!DOES NOT MANGAGE CODE RATES >0.74!####
###########################################
code_rate = 1 - 9/73
min = 0.1
max = 5
tol = 1e-7

while max - min > tol:
    sigma = (min + max)/2
    t, capacity = dmc_utils.optimize_thresholds(sigma, sigma, symmetric=True)
    result = capacity - code_rate
        
    if result > 0:
        min = sigma
    elif result < 0:
        max = sigma
    else:
        pass

sigma = (min + max)/2
rber = stats.norm.cdf(0, 1, sigma)
sigma_tot = sigma*2
print(f"""The shannon limit for the discretized A-AWGN channel ({code_rate:.4}, {0.5}) is {rber:.4}(RBER)/{sigma_tot:.4}(SIGMA)""")
#%% AWGN 
code_rate = 1 - 9/73
min = 0.1
max = 5
tol = 1e-7

while max - min > tol:
    sigma = (min + max)/2

    range = np.arange(-20,20,0.0001)
    cdf = 1/np.sqrt(8 * np.pi * sigma**2) * (np.exp(-(range-1)**2/(2*sigma**2)) + np.exp(-(range+1)**2/(2*sigma**2)))
    entropy = cdf * np.log2(np.sqrt(2*np.pi*np.e*sigma**2)*cdf)
    entropy[np.isnan(entropy)] = 0
    capacity = -np.trapz(entropy, range) 
    result = capacity - code_rate

    if result > 0:
        min = sigma
    elif result < 0:
        max = sigma
    else:
        raise Exception

sigma = (min + max)/2
rber = stats.norm.cdf(0, 1, sigma)
sigma_tot = sigma*2
print(f"The shannon limit for the continous AWGN channel ({code_rate:.4}) is {rber:.4}(RBER)/{sigma_tot:.4}(SIGMA)")

# %%
