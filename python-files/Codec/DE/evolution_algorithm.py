# %% Imports and variables
import time
import numpy as np
import density_evolution
import generate_distributions
import de_methods
import random

#%% Initial population generation
"""
One individual is represented by a n_vn * n_cn binary protograph matrix, from which the edge and node distributions may be calculated.
The initial generated by selecting 2 non zero indexes for each variable node, to be sure it has enough connections. 
"""
n_pop = 50
init_pop = np.random.randint()

#%% Fitness evaluation
"""
The inidvidual is evaluated first by a number of conditions:
    1. All rows and columns should sum to at least 2
    2. All rows should  sum to at most 20.

If the individual does not meet the conditions it is heavily punished, if it does it it is evaluated through density evolution.
"""

#%% Genetic methods
"""
In the creation of the next generation tournament selection is applied.

Three genetic processes are applied to the genome:
    1. Vertical crossover: Take a number of rows from each individual.
    2. Horizontal crossover: Similar to the vertical case but for rows
    3. Mutation: Flip a random tile of the matrix
"""

#%% Full algorithm