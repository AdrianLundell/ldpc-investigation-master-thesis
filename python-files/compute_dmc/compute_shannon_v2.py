import code
import dmc_utils
import yaml
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import math


def compute_pe(capacity,code_rate):
    if code_rate < capacity:
        return 0
    else:
        return 1-capacity/code_rate


def compute_capacity(sigma, skew):
    sigma1 = (1-skew)*sigma
    sigma2 = skew*sigma

    t, capacity = dmc_utils.optimize_thresholds(
        sigma1, sigma2, symmetric=True)
    return capacity

######


code_rate = 0.9
skew = 0.5
sigmas = np.arange(0.2, 1, 0.001)
pe = np.zeros(sigmas.size)
capacities = np.zeros(sigmas.size)
x = np.arange(-20, 20, 0.001)

for idx, sigma in enumerate(sigmas):
    capacity = compute_capacity(sigma,skew)
    capacities[idx] = capacity
    if capacity < code_rate:
        shannon_sigma = sigma
        break


cdf = stats.norm.cdf(0, 1, shannon_sigma)
rber = 100*cdf


##############
#Continuous
##############
code_rate = 0.87
sigmas = np.arange(0.25,1,0.001)
ebn0 = np.zeros(sigmas.size)
pe = np.zeros(sigmas.size)
capacities = np.zeros(sigmas.size)
x = np.arange(-20,20,0.0001)
for i,sigma in enumerate(sigmas):
    ebn0[i] = code_rate/(2*sigma**2)
    cdf = 1/np.sqrt(8*np.pi*sigma**2)*(np.exp(-(x-1)**2/(2*sigma**2))+np.exp(-(x+1)**2/(2*sigma**2)))
    
    entropy = cdf * np.log2(np.sqrt(2*np.pi*np.e*sigma**2)*cdf)
    entropy[np.isnan(entropy)] = 0
    capacity = -np.trapz(entropy, x)
    capacities[i] = capacity
    if capacity < code_rate:
        shannon_sigma = sigma
        break


cdf = stats.norm.cdf(0, 1, shannon_sigma)
rber = 100*cdf
print(rber)
#plt.plot(capacities, ebn0)
#plt.show()

