import dmc_utils
import yaml
import numpy as np 
import scipy.stats as stats 

skew = 0.5
code_rate = 1 - 9/73
is_symmetric = (skew == 0.5)

def compute_capacity(rber, skew):
    sigma = dmc_utils.rber_to_sigma(rber, skew)
    sigma1 = (1-skew)*sigma
    sigma2 = skew*sigma 

    t, capacity = dmc_utils.optimize_thresholds(sigma1, sigma2, symmetric = is_symmetric)
    return capacity

min = 0.001
max = 0.5
tol = 1e-7

while max - min > tol:
    x = (min + max)/2
    capacity = compute_capacity(x, skew) 
    result = (1-capacity) - code_rate

    if result < 0:
        min = x
    elif result > 0:
        max = x
    else:
        raise Exception

rber = (min + max)/2
sigma = dmc_utils.rber_to_sigma(rber, skew)
print(f"The shannon limit for the discretized A-AWGN channel ({code_rate}, {skew}) is {rber}(RBER)/{sigma}(SIGMA)")


min = 0.001
max = 0.5
tol = 1e-7

#Continous symmetric case
while max - min > tol:
    x = (min + max)/2

    s = dmc_utils.rber_to_sigma(x, skew)
    range = np.linspace(-10, 10, 1000)
    cdf = 1/np.sqrt(8 * np.pi * s**2) * (np.exp(-(range-1)**2/(2*s**2)) + np.exp(-(range+1)**2/(2*s**2)))
    entropy = cdf * np.log2(cdf)
    capacity = -np.trapz(entropy, range) - 1/2*np.log2(2*np.pi*np.e*s**2)
    result = (1-capacity) - code_rate

    if result < 0:
        min = x
    elif result > 0:
        max = x
    else:
        raise Exception

rber = (min + max)/2
sigma = dmc_utils.rber_to_sigma(rber, skew)
print(f"The shannon limit for the continous AWGN channel ({code_rate}) is {rber}(RBER)/{sigma}(SIGMA)")
