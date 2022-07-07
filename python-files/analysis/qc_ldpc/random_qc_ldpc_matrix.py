import os
import datetime
import numpy as np


def save(G,n,m,N, filename):
    """
    Saves the protograph as a .qc file, as defined here: https://aff3ct.readthedocs.io/en/latest/user/simulation/parameters/codec/ldpc/decoder.html
    Appends metadata at end
    """
    data = f"{n} {m} {N}\n\n"
    for row in G:
        for char in row:
            data = data + str(char) + " "
        data += "\n"
    
    with open(filename, "w") as f:
        f.write(data)
    with open(filename, "a") as f:
        metadata = f"""
# Random QC-code that is guaranteed to be invertible.
# 
# Name   : {filename}
# Metric : Minimum distance
# Date   : {datetime.datetime.now()}
"""
        f.write(metadata)


n = 73
m = 9
N = 256 # Scaling factor
output_file = "data/random.qc"

# Create generator matrix G that is guaranteed to be invertible
G = np.random.randint(2, size=(m, n))-1
I = np.identity(m,dtype=int)-1
G[-m:,-m:] = I

save(G, n, m, N, output_file)


