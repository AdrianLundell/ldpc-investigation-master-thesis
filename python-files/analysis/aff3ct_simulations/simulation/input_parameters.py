import numpy as np
from scipy import special


###################################################
# Convert interval of RBER to interval of Eb/N0
###################################################
rber_start = 0.0001
rber_end = 0.1
code_rate = 0.87671232876

def rber_to_sigma(rber): 
    sigma = -1 / (special.erfinv(2*rber - 1) * np.sqrt(2))
    return sigma

def sigma_to_ebn0(sigma):
    Es = 1
    n0 = 2*sigma**2
    snr = 1/sigma**2
    esn0 = 10 * np.log10(Es/n0)
    ebn0 = 10 * np.log10(Es/(code_rate*n0))
    return ebn0

sigma_start = rber_to_sigma(rber_start)
sigma_end = rber_to_sigma(rber_end)

ebn0_end = sigma_to_ebn0(sigma_start)
ebn0_start = sigma_to_ebn0(sigma_end)


message = f"""
RBER interval: {rber_start}-{rber_end}
Sigma interval: {sigma_start:.2f}-{sigma_end:.2f}
Eb/N0 interval (input into Aff3ct): {ebn0_start:.2f}-{ebn0_end:.2f}
"""

print(message)