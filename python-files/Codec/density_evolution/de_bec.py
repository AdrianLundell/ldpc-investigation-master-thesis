#%%
import numpy as np 
import matplotlib.pyplot as plt 
#%%

a = np.zeros(1024)
a[512-100:] = 0.1
a[512-50:] = 0.2
a[512+50:] = 0.5
a[512+100:] = 1
a = np.pad(a, (1024 + 512,1024 + 512), constant_values = (0,1))
A = np.fft.fft(a)

plt.plot(A)



#%% Binary erasure channel erasure limit
rho_coeffs = np.array([0, 0, 0, 0, 0, 1])
lambda_coeffs = np.array([0, 0, 1])
n_iter = 40

x0 = 1
xl = x0
dx = 0.0001

while xl > 10**-8:
    x0 = x0-dx
    xl = x0

    for i in range(n_iter):
        x1 = 1 - xl
        x_pow = np.array([x1**i for i in range(rho_coeffs.size)])
        x2 = sum(rho_coeffs * x_pow)
        x3 = 1 - x2
        x_pow = np.array([x3**i for i in range(lambda_coeffs.size)])
        xl = x0 * sum(lambda_coeffs * x_pow) 

print(x0)



#%% Binary erasure channel single run
rho_coeffs = np.array([0, 0, 0, 0, 0, 1])
lambda_coeffs = np.array([0, 0, 1])
n_iter = 10

x0 = 0.4
xl = x0

fig, axes = plt.subplots(2,2)
axes = np.reshape(axes, 4)

for i in range(n_iter):
    x1 = 1 - xl
    x_pow = np.array([x1**i for i in range(rho_coeffs.size)])
    x2 = sum(rho_coeffs * x_pow)
    x3 = 1 - x2
    x_pow = np.array([x3**i for i in range(lambda_coeffs.size)])
    xl = x0 * sum(lambda_coeffs * x_pow) 

    axes[0].scatter(i, x1)
    axes[1].scatter(i, x2)
    axes[2].scatter(i, x3)
    axes[3].scatter(i, xl)

# %%
