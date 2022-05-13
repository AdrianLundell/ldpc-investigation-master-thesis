#%%
import numpy as np
import scipy 
import scipy.special as sp 
import scipy.stats as stats
import matplotlib.pyplot as plt
import pandas as pd

#%% RBER to EBN0
code_rate = 64/73
rber = np.linspace(0.001, 0.1)
sigma = -1 / (sp.erfinv(2*rber - 1) * np.sqrt(2))
Es = 1
n0 = 2*sigma**2
snr = 1/sigma**2
esn0 = 10 * np.log10(Es/n0)
ebn0 = 10 * np.log10(Es/(code_rate*n0))

# %% Capacity 
x = np.linspace(-10, 10, 1000)
capacities = []
for s in sigma:
    cdf = 1/np.sqrt(8 * np.pi * s**2) * (np.exp(-(x-1)**2/(2*s**2)) + np.exp(-(x+1)**2/(2*s**2)))
    entropy = cdf * np.log2(cdf)
    capacity = -np.trapz(entropy, x) - 1/2*np.log2(2*np.pi*np.e*s**2)
    capacities.append(capacity)

#%%
plt.figure(figsize=(15,5))
plt.subplot(1,3,1)
plt.plot(ebn0, capacities)

df = pd.read_csv("Single_soft_capacity.csv", sep=", ", header=0)
df = df.drop_duplicates(subset='capacity', keep="last")
plt.plot(df.ebn0, df.capacity)

plt.grid()

plt.subplot(1,3,2)
plt.plot(rber, capacities)
plt.grid()

plt.subplot(1,3,3)
plt.plot(rber, ebn0)
plt.grid()
plt.plot()
# %% Try
mu1 = -1 
mu2 = 1

def mid_point():
    """Calculate the point of intersection for two gaussian distributions"""
    if sigma1 == sigma2:
        return (mu1+mu2)/2

    else:
        a = sigma2**2 - sigma1**2
        b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
        c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
        return (-b + np.sqrt(b**2 - 4*a*c))/(2*a)

sigma = 1.18818294989389
ratio = 0
sigma1 = sigma * (1+ratio)
sigma2 = sigma * (1-ratio)

p1 = stats.norm(-1,sigma1)
p2 = stats.norm(1, sigma2)

x1 = p1.rvs(100000)
x2 = p2.rvs(100000)
m = mid_point()
print(m)
print((sum(x1 > m) + sum(x2 < m))/200000)


sigma = 1
ratio = 0.1
sigma1 = 1.5684402773483785
sigma2 = 0.9055394163349377

p1 = stats.norm(-1, sigma1)
p2 = stats.norm(1, sigma2)

x1 = p1.rvs(100000)
x2 = p2.rvs(100000)
#m = mid_point()
print(m)
print((sum(x1 > m) + sum(x2 < m))/200000)

# %% Fixed RBER
ratio = 0.5
mu1 = -1
mu2 = 1
rber = 0.2
true_sigma = -1 / (sp.erfinv(2*rber - 1) * np.sqrt(2))

sigma = np.linspace(1, true_sigma+1,200)
rber = np.full(sigma.shape, rber) 
sigma1 = (1+ratio) * sigma 
sigma2 = (1-ratio) * sigma

if ratio == 0:
    t = (mu1+mu2)/2
else:
    a = sigma2**2 - sigma1**2
    b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
    c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
    t = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
 
rhs = 1/4*(2 + sp.erf((t-1)/(sigma1*np.sqrt(2))) - sp.erf((t+1)/(sigma2*np.sqrt(2))))
rhs2 = 1/2*(1 + sp.erf((t-1))/(sigma*np.sqrt(2)))

plt.plot(sigma, rber)
plt.plot(sigma, rhs)
#plt.plot(sigma, rhs2)
plt.plot([true_sigma,true_sigma],[0.1, 0.3], '--')
plt.show()

diff = np.abs(rber - rhs)
min_index = np.argmin(diff)
print(sigma[min_index], sigma1[min_index], sigma2[min_index], true_sigma)

# %% calculate variance 
sigma = 2
dsigma = 0.2
sigma1 = sigma*dsigma
sigma2 = sigma*(1-dsigma) 
N = 1000000

x1 = stats.norm.rvs(size = N, loc =-1, scale=sigma1)
x2 = stats.norm.rvs(size = N, loc =1, scale=sigma2)
x = np.append(x1, x2)

np.var(x)
# %%
# %% Fixed RBER
ratio = 0.5
mu1 = -1
mu2 = 1
rber = 0.02
true_sigma = -1 / (sp.erfinv(2*rber - 1) * np.sqrt(2))
sigma = np.linspace(0.01, true_sigma+2,200)


result = []

for ratio in np.linspace(0,1, endpoint=False):
    sigma1 = (ratio) * sigma 
    sigma2 = (1-ratio) * sigma

    if ratio == 0:
        t = (mu1+mu2)/2
    else:
        a = sigma2**2 - sigma1**2
        b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
        c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
        t = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    
    rhs = 1/4*(2 + sp.erf((t-1)/(sigma1*np.sqrt(2))) - sp.erf((t+1)/(sigma2*np.sqrt(2))))
    
    diff = np.abs(rber - rhs)
    min_index = np.argmin(diff)
    plt.plot(sigma, rhs)
    plt.show()
    result.append(sigma[min_index])

plt.plot(np.linspace(0,1), result)
# %%
