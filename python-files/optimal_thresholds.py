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

#%% Help functions
sigma1 = 0.9
sigma2 = 0.1
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

def plot(t):
    """Plot two gaussian distributions with thresholds t"""
    x = np.linspace(norm.ppf(0.01, loc=mu1, scale=sigma1), norm.ppf(0.99, loc=mu2, scale=sigma2), 100)
    
    plt.figure()
    plt.plot(x, norm.pdf(x, loc=mu1, scale=sigma1), 'black')
    plt.plot(x, norm.pdf(x, loc=mu2, scale=sigma2), 'black')
    for threshold in t:
        plt.plot([threshold, threshold], [0,0.7], '--')

    plt.show()

def LLR(t):
    """Calculate LLRs for BI-GAWN channel with thresholds t"""
    p = (np.array([[norm.cdf(t[0], mu1, sigma1), norm.cdf(t[1], mu1, sigma1), norm.cdf(t[2], mu1, sigma1), 1],
            [norm.cdf(t[0], mu2, sigma2), norm.cdf(t[1], mu2, sigma2),norm.cdf(t[2], mu2, sigma2), 1]]))
    
    yx_prob = np.array([p[0, 0], p[0, 1] - p[0, 0], p[0,2] - p[0, 1], p[0, 3] - p[0, 2], 
                        p[1, 0], p[1, 1] - p[1, 0], p[1,2] - p[1, 1], p[1, 3] - p[1, 2]])
    
    return  np.log(yx_prob[0:4]/yx_prob[4:])

#%% Calculate sigmas from noise level in DB for a channel with Es=1
noise_db = 4
N0 = 2/10**(noise_db/10)
sigma_ratio = 0.5 #Distribution of noise between the distributions

mu1 = -1
mu2 = 1
sigma1 = np.sqrt(N0*sigma_ratio)
sigma2 = np.sqrt(N0*(1-sigma_ratio))
m = mid_point()

plot([m])

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
m = mid_point()
positive_offset = 0.2
negative_offset = 0.2
t = [m-negative_offset, m, m+positive_offset]

plot(t)
print(f"MI(X,Y) = {mutual_info(t)}")

#%%Exhausive search of optimal mutual info for 1D problem of symmetric threshold
offsets = np.arange(0, 1912.5e-3, 7.5e-3)
mi_result = np.zeros(len(offsets))

sigma1 = 1
sigma2 = 1
m = mid_point()

print("Computing...")
for i, q in enumerate(offsets):
    t = [m-q, m, m+q]
    mi_result[i] = mutual_info(t)

plt.figure()
plt.plot(offsets, mi_result)
plt.show()

max_mi = np.max(mi_result)
max_mi_index = np.argmax(mi_result)

plot([m-offsets[max_mi_index], m, m+offsets[max_mi_index]])
print(f"MAX MI(X,Y) = {max_mi} for offset {offsets[max_mi_index]}")

print(LLR([-offsets[max_mi_index], 0, offsets[max_mi_index]]))

#%%Exhausive search of optimal mutual info for asymmetric thresholds
def optimize_threhsolds():
    m = middle
    offsets = np.arange(7.5e-3, 1, 7.5e-3)
    mi_result = np.zeros((len(offsets), len(offsets)))

    for i, positive_offset in enumerate(offsets):
        for j, negative_offset in enumerate(offsets):
            t = [m-negative_offset, m, m+positive_offset]
            mi_result[i, j] = mutual_info(t)

    #plt.figure()
    #plt.imshow(mi_result)
    #plt.colorbar()
    #plt.show()

    max_mi = np.max(mi_result)
    max_mi_index = np.unravel_index(np.argmax(mi_result), mi_result.shape)
    pos_offset = offsets[max_mi_index[0]]
    neg_offset = offsets[max_mi_index[1]]

    return max_mi, pos_offset, neg_offset

# %% Calculate many thresholds and output file with columns
#noise_db, sigma_ratio, sigma1, sigma2, max_mi, mid_point, neg_offset, pos_offset, llr1, llr2, llr3, llr4
db_range = np.arange(3,8,0.5)
sigma_range = [0.2,0.3,0.4,0.5]

result = np.zeros((len(db_range)*len(sigma_range), 12))
i = 0

for noise_db in db_range:
    for sigma_ratio in sigma_range:
        N0 = 2/10**(noise_db/10)
        sigma1 = np.sqrt(N0*sigma_ratio)
        sigma2 = np.sqrt(N0*(1-sigma_ratio))
        middle = mid_point()

        print(f"Calculating thresholds for SNR={noise_db}dB, sigma ratio={sigma_ratio}.")

        max_mi, pos_offset, neg_offset = optimize_threhsolds()
        llr = LLR(np.array([middle-neg_offset, middle, middle+pos_offset]))

        result[i,:] = np.array([noise_db, sigma_ratio, sigma1, sigma2, max_mi, middle, neg_offset, pos_offset, llr[0], llr[1], llr[2], llr[3]])

        i+=1

np.savetxt(fname = "Thresholds.csv", 
            X = result, 
            delimiter=", ",
            fmt='%f',
            header ="noise_db, sigma_ratio, sigma1, sigma2, max_mi, middle, neg_offset, pos_offset, llr1, llr2, llr3, llr4")
# %% Optimise mean performance of assymmetric thresholds 
db_range = np.arange(3,8,0.5)
sigma_range = [0.2,0.3,0.4,0.5]

offsets = np.arange(7.5e-3, 1, 7.5e-3)
mi_result = np.zeros((len(offsets), len(offsets)))

for i, positive_offset in enumerate(offsets):
    print(i)
    for j, negative_offset in enumerate(offsets):
        
        for noise_db in db_range:
            for sigma_ratio in sigma_range:
                N0 = 2/10**(noise_db/10)
                sigma1 = np.sqrt(N0*sigma_ratio)
                sigma2 = np.sqrt(N0*(1-sigma_ratio))
                middle = mid_point()

                t = [m-negative_offset, m, m+positive_offset]
                mi_result[i, j] += mutual_info(t)

max_mi = np.max(mi_result)/(len(db_range)*len(sigma_range))
max_mi_index = np.unravel_index(np.argmax(mi_result), mi_result.shape)

print(f"Max mean mutual info={max_mi} for offsets +{offsets[max_mi_index[0]]}, -{offsets[max_mi_index[1]]}")
# %% Save result
#Max mean mutual info=0.8432954142791612 for offsets +0.3825, -0.5175
# Optimise mean performance of symmetric thresholds 

db_range = np.arange(3,8,0.5)
sigma_range = [0.2,0.3,0.4,0.5]

offsets = np.arange(7.5e-3, 1, 7.5e-3)
mi_result = np.zeros(len(offsets))

for i, offset in enumerate(offsets):
        
    for noise_db in db_range:
        for sigma_ratio in sigma_range:
            N0 = 2/10**(noise_db/10)
            sigma1 = np.sqrt(N0*sigma_ratio)
            sigma2 = np.sqrt(N0*(1-sigma_ratio))
            middle = mid_point()

            t = [m-offset, m, m+offset]
            mi_result[i] += mutual_info(t)

max_mi = np.max(mi_result)/(len(db_range)*len(sigma_range))
max_mi_index = np.argmax(mi_result)

print(f"Max mean mutual info={max_mi} for offsets +-{offsets[max_mi_index]}")
# %%Save resuts
#Max mean mutual info=0.842982434826849 for offsets +-0.5099999999999999

#%% Compare different mus
snr_range = np.arange(4)
ratio_range = [0.4]
mu1_range = [-0.2, -0.5, -0.5, -1] 
mu2_range = [0.2, 0.5, 1, 2]
result = np.zeros((len(snr_range)*len(ratio_range)*len(mu1_range)*len(mu2_range), 12))
i = 0

for mu1 in mu1_range:
    for mu2 in mu2_range:
        for snr in snr_range:
            for ratio in ratio_range:
                
                sigma = np.sqrt((mu1-mu2)**2/10**(snr / 10))
                sigma1 = np.sqrt(ratio * sigma**2)
                sigma2 = np.sqrt((1-ratio) * sigma**2)
                
                middle = mid_point()

                print(f"Calculating thresholds for SNR={noise_db}dB, sigma ratio={sigma_ratio}.")

                max_mi, pos_offset, neg_offset = optimize_threhsolds()
                llr = LLR(np.array([middle-neg_offset, middle, middle+pos_offset]))

                result[i,:] = np.array([mu1, mu2, snr, ratio, max_mi, middle, neg_offset, pos_offset, llr[0], llr[1], llr[2], llr[3]])

                i+=1

np.savetxt(fname = "Thresholds_mu.csv", 
            X = result, 
            delimiter=", ",
            fmt='%f',
            header ="mu1, mu2, snr, ratio, max_mi, middle, neg_offset, pos_offset, llr1, llr2, llr3, llr4")
