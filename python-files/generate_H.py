#%% 
import numpy as np 
import scipy.linalg as sp 

C = 130
R = 69 


H = -1 * np.ones((R, C))

H[:,:R] = np.identity(R) - 1
H[:,R:C:10] = 1


for row in H:
    s = ''
    for nbr in row:
        s += str(int(nbr)) + " "
    print(s)

print(sp.lu(H))

# %%
