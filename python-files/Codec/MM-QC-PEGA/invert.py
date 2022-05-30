#%%
import numpy as np 
import galois 

proto = np.array([[ 0, -1,  1, -1, -1, -1,  4, -1, -1,  1, -1,  0, -1,  0, -1,  7, -1, -1,  8, -1],
       [-1,  7, -1, -1,  8,  5, -1, -1, -1,  0, -1,  9, -1,  8, -1, -1, -1, -1,  2,  6],
       [-1,  7, -1,  1, -1,  4, -1, -1,  4, -1,  8, -1,  5, -1, -1, -1,   5, -1, -1,  4],
       [-1, -1,  3, -1,  2, -1, -1,  7,  4, -1,  1, -1, -1, -1,  4, -1,    2,  1, -1, -1],
       [ 6, -1, -1,  9, -1, -1,  7,  1, -1, -1, -1, -1,  8, -1,  8,  0,  -1,  4, -1, -1]])
    
A = np.zeros(proto.shape)
A[proto > -1] = 1
A = galois.GF2(A.astype(int))
# %%
m = 5
n = 10
order = [0]
B = A[:,0:1]

for i in range(1,m): 
    for j in range(n):
        if j not in order:
            C = np.append(B, A[:,j:j+1], axis = -1)

            if np.linalg.matrix_rank(C) == i+1:
                B = C
                order.append(j)
                break

print("shape ", B.shape )         
print("rank ", np.linalg.matrix_rank(B))
print("det ", np.linalg.det(B))
print(order)
print("inv ", np.linalg.inv(B))
# %%
