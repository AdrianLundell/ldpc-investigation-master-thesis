#%% 
import numpy as np 

C = 581
R = 69 

H = np.zeros((R, C)) -1
H[:,:R] = np.identity(R) - 1
#H[:,R:] = np.random.randint(-1, 32-1, H[:,R:].shape)
H[:,R:2*R] =  np.identity(R) - 1
H[:,2*R:3*R] =  np.identity(R) - 1
H[:,3*R:4*R] =  np.identity(R) - 1
H[:,4*R:5*R] =  np.identity(R) - 1
H[:,5*R:6*R] =  np.identity(R) - 1
H[:,6*R:7*R] =  np.identity(R) - 1
H[:,7*R:8*R] =  np.identity(R) - 1
H[0:29,8*R:] =  np.identity(29) - 1



for row in H:
    s = ''
    for nbr in row:
        s += str(int(nbr)) + " "
    print(s)
