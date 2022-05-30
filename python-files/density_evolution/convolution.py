#%%
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as sp 

N = 256
n = 64
x1 = np.zeros(N)
x1[n:] = 1
x2 = np.zeros(N)
x2[n] = 1

y2 = sp.convolve(x1, x2)
plt.plot(y2)

#%%
y = sp.fftconvolve(x,x) * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)  * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)  * L/N
plt.plot(y)

N = 2048
x = np.zeros(N)
x[512:] = 1

y = sp.fftconvolve(x,x) * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)  * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)  * L/N
plt.plot(y)

# %%
#%%
N = 1024
x = np.zeros(N)
x[512:] = 1

y = sp.fftconvolve(x,x) * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)  * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)  * L/N
plt.plot(y)

plt.show()
N = 1024
x = np.zeros(N)
x[256:] = 1


y = sp.fftconvolve(x,x) * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)  * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)  * L/N
plt.plot(y)

# %%
N = 2048
x = np.zeros(N)
x[512 + 256:] = 1

y = sp.fftconvolve(x,x)[512:1024+512] * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)[512:1024+512]  * L/N
plt.plot(y)

y = sp.fftconvolve(x,y)[512:1024+512]  * L/N
plt.plot(y)

# %%
