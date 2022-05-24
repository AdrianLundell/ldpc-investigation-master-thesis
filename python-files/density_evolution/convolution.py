#%%
import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as sp 

N = 2048
X = np.arange(N)
L = 1
x1 = np.zeros(N)
x1[512:] = 1
x2 = np.zeros(N)
x2[512] = 1

X1 = np.fft.fft(x1)
X2 = np.fft.fft(x2)
y1 = np.abs(np.fft.ifft(X1 * X2))
y2 = sp.convolve(x1, x2)
#plt.plot(X1)
#plt.plot(X2)
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
