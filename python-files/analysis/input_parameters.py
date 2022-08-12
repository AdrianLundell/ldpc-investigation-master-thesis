#%%
import numpy as np
from scipy import special
from scipy.stats import norm

###################################################
# Convert interval of RBER to interval of Eb/N0
###################################################
rber_start = 0.0005
rber_end = 0.1
code_rate = 1 - 9/73

def rber_to_sigma(rber):
    sigma = -1 / (special.erfinv(2*rber - 1) * np.sqrt(2))
    return sigma

def test_sigma(sigma,rber):
    assert np.abs(norm.cdf(-1,loc=0,scale=sigma)-rber) < 1e-15, "Computed sigma and rber do not match"

def sigma_to_ebn0(sigma):
    Es = 1
    n0 = 2*sigma**2
    snr = 1/sigma**2
    esn0 = 10 * np.log10(Es/n0)
    ebn0 = 10 * np.log10(Es/(code_rate*n0))
    return ebn0

min_rber = 0.0005
max_rber = 0.1
n_simulations = 20

logged_range = np.linspace(np.log(min_rber), np.log(max_rber), n_simulations)

rber_range = np.exp(logged_range)
sigma_range = rber_to_sigma(rber_range)
ebn0_range = sigma_to_ebn0(sigma_range)

print(f"""
code rate = {code_rate}
RBER range : {repr(rber_range)}
EBN0 range : {repr(ebn0_range)}
""")

# sigma_start = rber_to_sigma(rber_start)
# sigma_end = rber_to_sigma(rber_end)

# ebn0_end, esn0_end = sigma_to_ebn0(sigma_start)
# ebn0_start, esn0_start = sigma_to_ebn0(sigma_end)

# # Test the outputed sigma
# test_sigma(sigma_start,rber_start)
# test_sigma(sigma_end,rber_end)

# for rber in np.linspace(rber_start, rber_end):
#     sigma = rber_to_sigma(rber)
#     ebn0 = sigma_to_ebn0(sigma)

#     print(f"{rber} {ebn0}")


# message = f"""
# RBER to EB/N0
# RBER interval: {rber_start}-{rber_end}
# Sigma interval: {sigma_start:.3f}-{sigma_end:.3f}
# Es/N0 interval: {esn0_start:.3f}-{esn0_end:.3f}
# Eb/N0 interval (input into Aff3ct): {ebn0_start:.3f}-{ebn0_end:.3f}
# """
# print(message)

# # Functions implemented as in Aff3ct
# def ebn0_to_esn0(ebn0,bit_rate,bps=1):
#     esn0 = ebn0+10*np.log10(bit_rate*bps)
#     return esn0

# def esn0_to_sigma(esn0,upsample_factor=1):
#     sigma = np.sqrt(upsample_factor/(2*10**(esn0/10)))
#     return sigma

# esn0_start = ebn0_to_esn0(ebn0_start,code_rate)
# esn0_end = ebn0_to_esn0(ebn0_end,code_rate)

# sigma_start = esn0_to_sigma(esn0_end)
# sigma_end = esn0_to_sigma(esn0_start)

# message = f"""
# Aff3ct Eb/N0 to sigma:
# Eb/N0 interval (input into Aff3ct): {ebn0_start:.3f}-{ebn0_end:.3f}
# Es/N0 interval: {esn0_start:.3f}-{esn0_end:.3f}
# Sigma interval: {sigma_start:.3f}-{sigma_end:.3f}
# """
# print(message)

# %% Output exact 