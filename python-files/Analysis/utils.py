#%%
"""Utility methods for result analysis"""
import numpy as np
import scipy.special as sp 
import scipy.stats as stats 

def rber_to_sigma(rbers, skew = 0.5, n_iter = 10, mu1 = -1, mu2 = 1):
    """
    Produces the map from raw bit error rates (RBER) to sigmas of the discretized (asymmetric) BIAWGN channel.
    RBER is here interpretted as the value of the overlapping tails of the two PDFs divided by two, corresponding to the 
    average probability of a faulty read from hard decoding. Sigmas are found as the roots to
        RBER - 1/2[(1-cdf1(t)) + cdf2(t)]
    found with Newtons method.

    Skew is the parameter defining the asymmetry if the channel. 
    """
    

    result = []
    for rber in rbers:
        for i in range(n_iter):
            
            #Good initial approximation
            sigma = skew*(4+np.log10(rber)) 
            
            #Calculate sigmas for each distribution
            sigma1 = (skew) * sigma 
            d_sigma1 = skew 
            sigma2 = (1-skew) * sigma
            d_sigma2 = (1-skew)

            if sigma1 == sigma2:
                #Return direct result for symmetric distributions
                sigma = -1 / (sp.erfinv(2*rber - 1) * np.sqrt(2)) *2
                break
            
            else:
                #Calculate middle point of distributions
                a = sigma2**2 - sigma1**2
                d_a = 2*sigma2*d_sigma2 - 2*sigma1*d_sigma1
                b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
                d_b = -(4*sigma2*d_sigma2*mu1 - 4*sigma1*d_sigma1*mu2)
                c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
                d_c = 2*sigma2*d_sigma2*mu1**2 - 2*sigma1*d_sigma1*mu2**2 - 4*sigma**3*skew*(1-skew)*np.log(sigma2**2/sigma1**2)
        
                e = b**2 - 4*a*c
                d_e = 2*b*d_b - 4*d_a*c - 4*a*d_c
                f = 2*a
                d_f = 2*d_a
                t = (-b + np.sqrt(e))/(f)
                d_t = -d_b + (0.5/(np.sqrt(e))*d_e*f - np.sqrt(e)*d_f)/f**2
        
                #Calculate arguments to erf-funcions
                z1 = (t-1)/(sigma1*np.sqrt(2))
                z2 = (t+1)/(sigma2*np.sqrt(2))
                d_z1 = (d_t*sigma1 - (t-1)*d_sigma1)/(sigma1**2*np.sqrt(2))
                d_z2 = (d_t*sigma2 - (t+1)*d_sigma2)/(sigma2**2*np.sqrt(2))
                
                #Compute function value and derivative
                y = rber - 1/4*(2 + sp.erf(z1) - sp.erf(z2))
                d_y = -1/4*(d_z1 * 2/np.sqrt(np.pi) * np.exp(-z1**2) - d_z2 * 2/np.sqrt(np.pi) * np.exp(-z2**2))
                
                #Newtons method
                sigma = sigma - y/d_y

        result.append(sigma)

    return np.array(result)

#%%
def mid_point():
    """Calculate the point of intersection for two gaussian distributions"""
    mu1 = -1
    mu2 = 1
    if sigma1 == sigma2:
        return (mu1+mu2)/2

    else:
        a = sigma2**2 - sigma1**2
        b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
        c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
        return (-b + np.sqrt(b**2 - 4*a*c))/(2*a)

skew = 0.4
sigma = rber_to_sigma([0.01], skew)
sigma1 = sigma * (skew)
sigma2 = sigma * (1-skew)

p1 = stats.norm(-1,sigma1)
p2 = stats.norm(1, sigma2)
x1 = p1.rvs(100000)
x2 = p2.rvs(100000)
m = mid_point()
print((sum(x1 > m) + sum(x2 < m))/200000)

# %%
