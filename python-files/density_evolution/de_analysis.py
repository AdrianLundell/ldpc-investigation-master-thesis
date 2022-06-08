
# %% Imports
import numpy as np
import matplotlib.pyplot as plt
import math
from pyrsistent import v
from sklearn.feature_selection import f_classif
from numba import types
from numba import jit
from scipy.special import binom
from time import process_time
from density_evolution_symmetric import de_algorithm


# %%
# Degree distributions
lam = np.array([0, 0.2895, 0.3158, 0, 0, 0.3947])
rho = np.array([0, 0, 0, 0, 0, 0.9032, 0.0968])

# Sample channel (probability density of messages)
# Represents ext = np.arange(-30,30.01, 0.01)

# Extrinsic messages generator [min,step,size]
ext_min = -30
ext_step = 0.01
ext_size = 6001
p_0 = np.zeros(ext_size)
p_0[3239] = 91.606
p_0[2761] = 8.394
# Notice that the integral over the channel becomes 1

# Extrinsic messages are mapped onto gamma [min,step,size]
gamma_min = -10
gamma_step = 0.0002
gamma_size = 50000

# Density evolution
max_iter = 50

alg = de_algorithm(p_0, ext_min, ext_step, ext_size, gamma_min,
                   gamma_step, gamma_size, lam, rho, max_iter)
alg.density_evolution()

print(alg.p_e)
print(process_time())
