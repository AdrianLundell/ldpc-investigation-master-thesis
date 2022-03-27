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
#Not necessary, just a matter of rescaling and moving.
def mid_point(mu1, sigma1, mu2, sigma2):
    if sigma1 == sigma2:
        return (mu1+mu2)/2

    else:
        a = sigma2**2 - sigma1**2
        b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
        c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
        return (-b + np.sqrt(b**2 - 4*a*c))/(2*a)

def plot(t):
    x = np.linspace(norm.ppf(0.01, loc=-1, scale=sigma1), norm.ppf(0.99, loc=1, scale=sigma2), 100)
    
    plt.figure()
    plt.plot(x, norm.pdf(x, loc=-1, scale=sigma1), 'black')
    plt.plot(x, norm.pdf(x, loc=1, scale=sigma2), 'black')
    for threshold in t:
        plt.plot([threshold, threshold], [0,0.7], '--')

    plt.show()

#%% Calculate sigmas from noise level in DB for a channel with Es=1
noise_db = 4
N0 = 2/10**(noise_db/10)
sigma_ratio = 0.5 #Distribution of noise between the distributions

sigma1 = np.sqrt(N0*sigma_ratio)
sigma2 = np.sqrt(N0*(1-sigma_ratio))

plot([0])

#%% Calculating the mutual information between two gaussian distributions
def mutual_info(t):   
    p = (np.array([[norm.cdf(t[0], -1, sigma1), norm.cdf(t[1], -1, sigma1), norm.cdf(t[2], -1, sigma1), 1],
            [norm.cdf(t[0], 1, sigma2), norm.cdf(t[1], 1, sigma2),norm.cdf(t[2], 1, sigma2), 1]]))
    
    yx_prob = np.array([p[0, 0], p[0, 1] - p[0, 0], p[0,2] - p[0, 1], p[0, 3] - p[0, 2], 
                        p[1, 0], p[1, 1] - p[1, 0], p[1,2] - p[1, 1], p[1, 3] - p[1, 2]])

    y_prob = np.array([(yx_prob[0]+yx_prob[4])*0.5, (yx_prob[1]+yx_prob[5])*0.5, 
                       (yx_prob[2]+yx_prob[6])*0.5, (yx_prob[3]+yx_prob[7])*0.5])

    yx_prob = yx_prob[yx_prob>0]
    y_prob = y_prob[y_prob>0]

    y_entropy = -(y_prob @ np.log2(y_prob))
    yx_entropy = -(yx_prob @ np.log2(yx_prob))

    #Note: H(Y|X) = 1/2H(Y|X=0) + 1/2H(Y|X=1)
    #here both terms are concatenated.
    return y_entropy - yx_entropy*0.5

#Demonstration
positive_offset = 0.2
negative_offset = -0.2
t = [negative_offset, 0, positive_offset]

plot(t)
print(f"MI(X,Y) = {mutual_info(t)}")

#%%Exhausive search of optimal mutual info for 1D problem of symmetric threshold
offsets = np.arange(0, 1912.5e-3, 7.5e-3)
mi_result = np.zeros(len(offsets))

print("Computing...")
for i, q in enumerate(offsets):
    t = [-q, 0, q]
    mi_result[i] = mutual_info(t)

plt.figure()
plt.plot(offsets, mi_result)
plt.show()

max_mi = np.max(mi_result)
max_mi_index = np.argmax(mi_result)

plot([-offsets[max_mi_index], 0, offsets[max_mi_index]])
print(f"MAX MI(X,Y) = {max_mi} for offset {offsets[max_mi_index]}")


#%%Exhausive search of optimal mutual info for asymmetric thresholds
offsets = np.arange(7.5e-3, 1, 7.5e-3)
mi_result = np.zeros((len(offsets), len(offsets)))

print("Computing...")
for i, positive_offset in enumerate(offsets):
    print(f"{i}/{len(offsets)}")
    for j, negative_offset in enumerate(offsets):
        t = [0-negative_offset, 0, 0+positive_offset]
        mi_result[i, j] = mutual_info(t)

plt.figure()
plt.imshow(mi_result)
plt.colorbar()
plt.show()

max_mi = np.max(mi_result)
max_mi_index = np.unravel_index(np.argmax(mi_result), mi_result.shape)

plot([-offsets[max_mi_index[0]], 0, offsets[max_mi_index[1]]])
print(f"MAX MI(X,Y) = {max_mi} for offsets {offsets[max_mi_index[0]]}, {offsets[max_mi_index[1]]}")
