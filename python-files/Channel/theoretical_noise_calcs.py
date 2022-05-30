#%%
import re
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

#%%
def plot(t, sigma1, sigma2, mu1 = -1, mu2 =1):
    """Plot two gaussian distributions with thresholds t"""
    x = np.linspace(stats.norm.ppf(0.01, loc=mu1, scale=sigma1), stats.norm.ppf(0.99, loc=mu2, scale=sigma2), 100)
    
    plt.figure()
    plt.plot(x, stats.norm.pdf(x, loc=mu1, scale=sigma1), 'black')
    plt.plot(x, stats.norm.pdf(x, loc=mu2, scale=sigma2), 'black')
    for threshold in t:
        plt.plot([threshold, threshold], [0,0.7], '--')

    plt.show()
# %% Fixed RBER
ratio = 0.4
mu1 = -1
mu2 = 1
rber = 0.01
true_sigma = -1 / (sp.erfinv(2*rber - 1) * np.sqrt(2))*2

sigma = np.linspace(0.1, true_sigma+5,200)
rber = np.full(sigma.shape, rber) 
sigma1 = np.sqrt(ratio) * sigma 
sigma2 = np.sqrt(1-ratio) * sigma

if ratio == 0.5:
    t = np.array([(mu1+mu2)/2]*200)
else:
    a = sigma2**2 - sigma1**2
    b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
    c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
    t = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
 
rhs = 1/4*(2 + sp.erf((t-1)/(sigma1*np.sqrt(2))) - sp.erf((t+1)/(sigma2*np.sqrt(2))))

plt.plot(sigma, rber)
plt.plot(sigma, rhs)
plt.plot([true_sigma,true_sigma],[0.1, 0.3], '--')
plt.show()

diff = np.abs(rber - rhs)
min_index = np.argmin(diff)
print(f"sigma, sigma1, sigma2: ", sigma[min_index], sigma1[min_index], sigma2[min_index])
print(f"Theoretical symmetric: {true_sigma}")
plot([t[min_index]], sigma1[min_index],sigma2[min_index])

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
mu1 = -1
mu2 = 1
rber = 0.01
sigma = np.linspace(0.01, true_sigma+10,500)

result = []
temp = []
ratios = np.linspace(0.2,0.5, endpoint=True, num = 10)

for ratio in ratios:
    sigma1 = (1-ratio) * sigma 
    sigma2 = (ratio) * sigma

    if ratio == 0.5:    
        t = 0    

    else:
        a = sigma2**2 - sigma1**2
        b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
        c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
        t = (-b + np.sqrt(b**2 - 4*a*c))/(2*a)
    
    rhs = 1/4*(2 + sp.erf((t-mu2)/(sigma1*np.sqrt(2))) - sp.erf((t-mu1)/(sigma2*np.sqrt(2))))
    
    diff = np.abs(rber - rhs)
    min_index = np.argmin(diff)
    #plt.plot(sigma, rhs)
    #plt.show()
    result.append(sigma[min_index])

    if not ratio==0.5:
        t = t[min_index]
    sigma1 = sigma1[min_index]
    sigma2 = sigma2[min_index]

    p1 = 1 - stats.norm.cdf(t, loc=-1, scale = sigma1)
    p2 = stats.norm.cdf(t, loc=1, scale = sigma2)
    temp.append(p1 + p2)    

plt.plot(ratios, temp)
plt.show()
plt.plot(ratios, result)
# %%
# %% Fixed RBER
ratio = 0.5
mu1 = -1
mu2 = 1
rber = 0.01
#sigma = np.linspace(0.01, true_sigma+2,200)
N = 10
ratios = np.linspace(0.2,0.5, endpoint=True, num=N)
result = []
for ratio in ratios:

    sigma = ratio*(4+np.log10(rber))

    for i in range(10):
        sigma1 = (ratio) * sigma 
        d_sigma1 = ratio 
        sigma2 = (1-ratio) * sigma
        d_sigma2 = (1-ratio)

        if sigma1 == sigma2:
            sigma = -1 / (sp.erfinv(2*rber - 1) * np.sqrt(2)) * 2
            break
        
        else:
            a = sigma2**2 - sigma1**2
            d_a = 2*sigma2*d_sigma2 - 2*sigma1*d_sigma1
            b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
            d_b = -(4*sigma2*d_sigma2*mu1 - 4*sigma1*d_sigma1*mu2)
            c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
            d_c = 2*sigma2*d_sigma2*mu1**2 - 2*sigma1*d_sigma1*mu2**2 - 4*sigma**3*ratio*(1-ratio)*np.log(sigma2**2/sigma1**2)
    
            e = b**2 - 4*a*c
            d_e = 2*b*d_b - 4*d_a*c - 4*a*d_c
            f = 2*a
            d_f = 2*d_a
            t = (-b + np.sqrt(e))/(f)
            d_t = -d_b + (0.5/(np.sqrt(e))*d_e*f - np.sqrt(e)*d_f)/f**2
    
            z1 = (t-1)/(sigma1*np.sqrt(2))
            z2 = (t+1)/(sigma2*np.sqrt(2))
            d_z1 = (d_t*sigma1 - (t-1)*d_sigma1)/(sigma1**2*np.sqrt(2))
            d_z2 = (d_t*sigma2 - (t+1)*d_sigma2)/(sigma2**2*np.sqrt(2))
            
            y = 1/4*(2 + sp.erf(z1) - sp.erf(z2)) - rber
            d_y = 1/4*(d_z1 * 2/np.sqrt(np.pi) * np.exp(-z1**2) - d_z2 * 2/np.sqrt(np.pi) * np.exp(-z2**2))
            sigma = sigma - y/d_y

    result.append(sigma)

plt.plot(ratios,result)
# %%
