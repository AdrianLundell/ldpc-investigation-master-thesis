"""
This script implements optimal threshold setting with known voltage distributions for the asymmetric BI-AWGN distributions by maximising
the mutual information over the resulting DMC channel for stochastic variables input X and output Y depending on T 
    I(X,Y(T)) = H(Y(T)) - H(Y(T)|X)
With H denoting the entropy function.
"""

#%%
from scipy.stats import norm
import numpy as np 
import matplotlib.pyplot as plt 

#%% Calculating the optimal midpoint between two gaussian distributions
def mid_point(mu1, sigma1, mu2, sigma2):
    if sigma1 == sigma2:
        return (mu1+mu2)/2

    else:
        a = sigma2**2 - sigma1**2
        b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
        c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
        return (-b + np.sqrt(b**2 - 4*a*c))/(2*a)

#Demonstration
mu1 = -1
sigma1 = 0.25
mu2 =  1
sigma2 = 0.25
m = mid_point(mu1, sigma1, mu2, sigma2)

plt.figure()
x = np.linspace(norm.ppf(0.01, loc=mu1, scale=sigma1), norm.ppf(0.99, loc=mu2, scale=sigma2), 100)
plt.plot(x, norm.pdf(x, loc=mu1, scale=sigma1), 'black')
plt.plot(x, norm.pdf(x, loc=mu2, scale=sigma2), 'black')
plt.plot([m, m], [0,0.7], '--')

#%% Calculating the mutual information between two gaussian distributions
def mutual_info(t):
    p = np.exp(np.array([[norm.cdf(t[0], mu1, sigma1), norm.cdf(t[1], mu1, sigma1), norm.cdf(t[2], mu1, sigma1), 1],
            [norm.cdf(t[0], mu2, sigma2), norm.cdf(t[1], mu2, sigma2),norm.cdf(t[2], mu2, sigma2), 1]]))
    
    yx_prob = np.log(np.array([p[0, 0], p[0, 1]/p[0, 0], p[0,2]/p[0, 1], p[0, 3]/p[0, 2], 
                        p[1, 0], p[1, 1]/p[1, 0], p[1,2]/p[1, 1], p[1, 3]/p[1, 2]]))
    
    y_prob = np.log(np.array([(yx_prob[0]*yx_prob[4])**0.5, (yx_prob[1]*yx_prob[5])**0.5, 
                       (yx_prob[2]*yx_prob[6])**0.5, (yx_prob[3]*yx_prob[7])**0.5]))
    
    yx_prob = yx_prob[yx_prob==0]
    y_prob = y_prob[y_prob==0]

    y_entropy = -(np.sum(y_prob * np.log2(y_prob)))
    yx_entropy = -(np.sum(yx_prob * np.log2(yx_prob)))

    return y_entropy - yx_entropy

def no_log_mutual_info(t):
    p = (np.array([[norm.cdf(t[0], mu1, sigma1), norm.cdf(t[1], mu1, sigma1), norm.cdf(t[2], mu1, sigma1), 1],
            [norm.cdf(t[0], mu2, sigma2), norm.cdf(t[1], mu2, sigma2),norm.cdf(t[2], mu2, sigma2), 1]]))
    
    yx_prob = np.array([p[0, 0], p[0, 1] - p[0, 0], p[0,2] - p[0, 1], p[0, 3] - p[0, 2], 
                        p[1, 0], p[1, 1] - p[1, 0], p[1,2] - p[1, 1], p[1, 3] - p[1, 2]])
    
    y_prob = np.array([(yx_prob[0]*yx_prob[4])**0.5, (yx_prob[1]*yx_prob[5])**0.5, 
                       (yx_prob[2]*yx_prob[6])**0.5, (yx_prob[3]*yx_prob[7])**0.5])
    
    y_entropy = -(y_prob @ np.log2(y_prob))
    yx_entropy = -(yx_prob @ np.log2(yx_prob))

    return y_entropy - yx_entropy

#Demonstration
positive_offset = 0.2
negative_offset = 0.2
t = [m-negative_offset, m, m+positive_offset]

plt.figure()
x = np.linspace(norm.ppf(0.01, loc=mu1, scale=sigma1), norm.ppf(0.99, loc=mu2, scale=sigma2), 100)
plt.plot(x, norm.pdf(x, loc=mu1, scale=sigma1), 'black')
plt.plot(x, norm.pdf(x, loc=mu2, scale=sigma2), 'black')
plt.plot([m, m], [0,0.7], '--')
plt.plot([t[0], t[0]], [0,0.7], '--')
plt.plot([t[2], t[2]], [0,0.7], '--')
print(f"MI(X,Y) = {mutual_info(t)}")


#%%Exhausive search of optimal mutual info in relevant range
#offsets = np.arange(0, 1912.5e-3, 7.5e-3)
offsets = np.arange(7.5e-3, 1, 0.2)
mi_result = np.zeros((len(offsets), len(offsets)))

print("Computing...")
for i, positive_offset in enumerate(offsets):
    print(f"{i}/{len(offsets)}")
    for j, negative_offset in enumerate(offsets):
        t = [m-negative_offset, m, m+positive_offset]
        mi_result[i, j] = mutual_info(t)

plt.figure()
plt.imshow(mi_result)
plt.colorbar()
plt.show()

max_mi = np.min(mi_result)
max_mi_index = np.unravel_index(np.argmin(mi_result), mi_result.shape)

plt.figure()
x = np.linspace(norm.ppf(0.01, loc=mu1, scale=sigma1), norm.ppf(0.99, loc=mu2, scale=sigma2), 100)
plt.plot(x, norm.pdf(x, loc=mu1, scale=sigma1), 'black')
plt.plot(x, norm.pdf(x, loc=mu2, scale=sigma2), 'black')
plt.plot([m, m], [0,0.7], '--')
plt.plot([m + offsets[max_mi_index[0]], m + offsets[max_mi_index[0]]], [0,0.7], '--')
plt.plot([m - offsets[max_mi_index[1]], m - offsets[max_mi_index[1]]], [0,0.7], '--')

print(f"MAX MI(X,Y) = {max_mi} for offsets {offsets[max_mi_index[0]]}, {offsets[max_mi_index[1]]}")

#%%Exhausive search of optimal mutual info for symmetric threshold
#offsets = np.arange(0, 1912.5e-3, 7.5e-3)
offsets = np.arange(7.5e-3, 3, 0.2)
mi_result = np.zeros(len(offsets))

print("Computing...")
for i, positive_offset in enumerate(offsets):
    t = [m-positive_offset, m, m+positive_offset]
    mi_result[i] = mutual_info(t)

plt.figure()
plt.plot(offsets, mi_result)
plt.show()

max_mi = np.max(mi_result)
max_mi_index = np.argmax(mi_result)

plt.figure()
x = np.linspace(norm.ppf(0.01, loc=mu1, scale=sigma1), norm.ppf(0.99, loc=mu2, scale=sigma2), 100)
plt.plot(x, norm.pdf(x, loc=mu1, scale=sigma1), 'black')
plt.plot(x, norm.pdf(x, loc=mu2, scale=sigma2), 'black')
plt.plot([m, m], [0,0.7], '--')
plt.plot([m + offsets[max_mi_index], m + offsets[max_mi_index]], [0,0.7], '--')
plt.plot([m - offsets[max_mi_index], m - offsets[max_mi_index]], [0,0.7], '--')

print(f"MAX MI(X,Y) = {max_mi} for offset {offsets[max_mi_index]}")