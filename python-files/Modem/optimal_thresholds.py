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
sigma1 = 0.5
sigma2 = 0.5
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
    t = list(t)
    t.append(np.inf)
    p = np.array([norm.cdf(np.array(t), mu1, sigma1),
                  norm.cdf(np.array(t), mu2, sigma2)])
    
    yx_prob = np.array([np.ediff1d(p[0,:], to_begin = p[0,0]),
                        np.ediff1d(p[1,:], to_begin = p[1,0])])
    
    return  np.log(yx_prob[0,:]/yx_prob[1,:])

#%%
mu1 = -1 
mu2 = 1
sigma1 = 1
sigma2 = 1

x = mid_point()
print(1 - norm.cdf(x, loc=mu1, scale=sigma1))
print(norm.cdf(x, loc=mu2, scale=sigma2))


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

#%%Exhausive search of optimal mutual info
def optimize_threhsolds(symmetric = False, offsets = np.arange(7.5e-3, 1, 7.5e-3)):
    m = middle
    mi_result = np.full((len(offsets), len(offsets)), -np.inf)

    for i, positive_offset in enumerate(offsets):
        if symmetric:
            t = [m-positive_offset, m, m+positive_offset]
            mi_result[i, i] = mutual_info(t)
        else:
            for j, negative_offset in enumerate(offsets):
                t = [m-negative_offset, m, m+positive_offset]
                mi_result[i, j] = mutual_info(t)

    max_mi = np.max(mi_result)
    max_mi_index = np.unravel_index(np.argmax(mi_result), mi_result.shape)
    pos_offset = offsets[max_mi_index[0]]
    neg_offset = offsets[max_mi_index[1]]

    return max_mi, pos_offset, neg_offset

# %% Calculate many thresholds and output file with column
#ebn0, ratio, T0, llr1, llr2
ebn0_range = np.arange(0,10,0.5)
ratio_range = [0,1]
code_rate = 64/73
mu1 = -1
mu2 = 1

result = np.zeros((len(ebn0_range)*len(ratio_range), 6))
i = 0
for ebn0 in ebn0_range:
    for sigma_ratio in ratio_range:

        esn0 = ebn0 + 10*np.log10(code_rate)
        sigma = np.sqrt(1/(2*10**(esn0/10)))
        sigma1 = sigma
        sigma2 = sigma

        middle = mid_point()

        print(f"Calculating thresholds for ebn0={ebn0}dB, sigma ratio={sigma_ratio}.")
        t0 = middle
        llr = LLR(np.array([t0]))
        capacity = mutual_info(t0)
        result[i,:] = np.array([ebn0, sigma_ratio, t0, llr[0], llr[1], capacity])

        i+=1

np.savetxt(fname = "Thresholds_symmetric_hard.csv", 
            X = result, 
            delimiter=", ",
            fmt='%f',
            header ="ebn0, sigma_ratio, t0, llr0, llr1, capacity")

# %% Calculate many thresholds and output file with columns
#ebn0, ratio, T0, T1, T2, llr1, llr2, llr3, llr4
ebn0_range = np.arange(0,10,0.5)
ratio_range = [0,1]
code_rate = 64/73
mu1 = -1
mu2 = 1

result = np.zeros((len(ebn0_range)*len(ratio_range), 10))
i = 0
for ebn0 in ebn0_range:
    for sigma_ratio in ratio_range:

        esn0 = ebn0 + 10*np.log10(code_rate)
        sigma = np.sqrt(1/(2*10**(esn0/10)))
        sigma1 = sigma
        sigma2 = sigma

        middle = mid_point()

        print(f"Calculating thresholds for ebn0={ebn0}dB, sigma ratio={sigma_ratio}.")

        max_mi, pos_offset, neg_offset = optimize_threhsolds(symmetric=True)
        t0 = middle - neg_offset
        t1 = middle
        t2 = middle + neg_offset
        llr = LLR(np.array([t0, t1, t2]))
        c = mutual_info([t0,t1,t2])
        result[i,:] = np.array([ebn0, sigma_ratio, t0, t1, t2, llr[0], llr[1], llr[2], llr[3], c])

        i+=1

np.savetxt(fname = "Single_soft_capacity.csv", 
            X = result, 
            delimiter=", ",
            fmt='%f',
            header ="ebn0, sigma_ratio, t0, t1, t2, llr0, llr1, llr2, llr3, capacity")

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

db_range = np.arange(0,8,0.5)
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
