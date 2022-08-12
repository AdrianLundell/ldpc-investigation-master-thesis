#%%
import os
import datetime
import numpy as np
import galois 

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
# Date   : {datetime.datetime.now()}
"""
        f.write(metadata)


n_cn = 9
n_vn = 73
N = 256 # Scaling factor
output_file = "data/random.qc"

# Create generator matrix G that is guaranteed to be invertible
while True:
    proto = np.zeros((n_cn, n_vn), int)
    for j in range(n_vn):
        connections = np.random.randint(2,4)
        k = np.random.choice(n_cn, connections, replace=False)
        proto[k, j] = 1

    original = np.copy(proto)

    if np.linalg.matrix_rank(galois.GF2(proto)) == n_cn: 
        indexes = np.arange(proto.shape[1])

        for i in indexes:
            temp = proto.copy()
            temp[:, i] = 0
        
            if np.linalg.matrix_rank(temp) == proto.shape[0]:
                proto = temp
                indexes[i] = -1

        indexes = indexes[indexes >= 0]
        reminding_columns = [i for i in range(n_vn) if i not in indexes]
        new_order = reminding_columns + list(indexes)

        proto = np.take(original, new_order, 1)

        x = galois.GF2(proto[:,:n_cn])
        try:
            np.linalg.inv(x)
            proto = proto-1
            random = np.random.randint(257, size=(n_cn, n_vn))
            proto[proto==0] = random[proto == 0]

            save(proto, n_vn, n_cn, N, output_file)
            break
        except:
            pass



# %%
