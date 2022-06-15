import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import scipy.special as special

def mid_point(sigma1, sigma2, mu1 = -1, mu2 = 1):
    """Calculate the point of intersection for two gaussian distributions"""
    if sigma1 == sigma2:
        return (mu1+mu2)/2

    else:
        a = sigma2**2 - sigma1**2
        b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
        c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
        return (-b + np.sqrt(b**2 - 4*a*c))/(2*a)

def llrs(t, sigma1, sigma2, mu1 = -1, mu2 = 1):
    """Calculate LLRs for the A-BIAWGN channel with thresholds t"""
    t = list(t)
    t.append(np.inf)
    p = np.array([stats.norm.cdf(np.array(t), mu1, sigma1),
                  stats.norm.cdf(np.array(t), mu2, sigma2)])
    
    yx_prob = np.array([np.ediff1d(p[0,:], to_begin = p[0,0]),
                        np.ediff1d(p[1,:], to_begin = p[1,0])])
    
    return  np.log(yx_prob[0,:]/yx_prob[1,:])

def mutual_info(t, sigma1, sigma2, mu1 = -1, mu2 = 1):
    """Calculates the mutual infromation of the A-BIAWGN channel discretized with thresholds t"""
    p = np.array([stats.norm.cdf(np.array(t), mu1, sigma1),
                  stats.norm.cdf(np.array(t), mu2, sigma2)])
    
    yx_prob = np.array(([np.ediff1d(p[0,:], to_begin = p[0,0])],
                       [np.ediff1d(p[1,:], to_begin = p[1,0])]))

    y_prob = 0.5 * (yx_prob[0,:] + yx_prob[1,:])

    y_prob = y_prob.flatten()
    y_entropy = -(y_prob @ np.log2(y_prob))

    yx_prob = yx_prob.flatten()
    yx_entropy = -(yx_prob @ np.log2(yx_prob))*0.5

    return y_entropy - yx_entropy

def optimize_thresholds(sigma1, sigma2, mu1 = -1, mu2 = 1, symmetric = False, offsets = np.arange(7.5e-3, 1, 7.5e-3)):
    """Returns optimal offsets for three thresholds
    TODO: Generalize to arbitrary amount of thresholds
    TODO: More clever optimization strategy
    """
    m = mid_point(sigma1, sigma2, mu1, mu2)

    mi_result = np.full((len(offsets), len(offsets)), -np.inf)

    for i, positive_offset in enumerate(offsets):
        if symmetric:
            t = [m-positive_offset, m, m+positive_offset]
            mi_result[i, i] = mutual_info(t, sigma1, sigma2, mu1, mu2)
        else:
            for j, negative_offset in enumerate(offsets):
                t = [m-negative_offset, m, m+positive_offset]
                mi_result[i, j] = mutual_info(t, sigma1, sigma2, mu1, mu2)

    max_mi = np.max(mi_result)
    max_mi_index = np.unravel_index(np.argmax(mi_result), mi_result.shape)
    pos_offset = offsets[max_mi_index[0]]
    neg_offset = offsets[max_mi_index[1]]

    return np.array([m-neg_offset, m, m+pos_offset]), max_mi

def rber_to_sigma(rber, skew = 0.5, n_iter = 10, mu1 = -1, mu2 = 1):
    """
    Produces the map from raw bit error rates (RBER) to sigmas of the discretized (asymmetric) BIAWGN channel.
    RBER is here interpretted as the value of the overlapping tails of the two PDFs divided by two, corresponding to the 
    average probability of a faulty read from hard decoding. Sigmas are found as the roots to
        RBER - 1/2[(1-cdf1(t)) + cdf2(t)]
    using Newtons method.

    Skew is the parameter defining the asymmetry of the channel. 
    """

    assert 0.0001 <= rber <= 0.5, "RBER out of range."

    #Good initial approximation
    sigma = (1-skew)*(4+np.log10(rber)) 
    
    for i in range(n_iter):
        
        #Calculate sigmas for each distribution
        sigma1 = (1-skew) * sigma 
        d_sigma1 = (1-skew) 
        sigma2 = (skew) * sigma
        d_sigma2 = (skew)

        if skew == 0.5:
            #Return direct result for symmetric distributions
            sigma = -1 / (special.erfinv(2*rber - 1) * np.sqrt(2)) * 2
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
            z1 = (t-mu2)/(sigma2*np.sqrt(2))
            z2 = (t-mu1)/(sigma1*np.sqrt(2))
            d_z1 = (d_t*sigma2 - (t-mu2)*d_sigma2)/(sigma2**2*np.sqrt(2))
            d_z2 = (d_t*sigma1 - (t-mu1)*d_sigma1)/(sigma1**2*np.sqrt(2))
            
            #Compute function value and derivative
            y = 1/4*(2 + special.erf(z1) - special.erf(z2)) - rber
            d_y = 1/4*(d_z1 * 2/np.sqrt(np.pi) * np.exp(-z1**2) - d_z2 * 2/np.sqrt(np.pi) * np.exp(-z2**2))

            #Newtons method
            sigma = sigma - y/d_y

    return sigma

def plot_awgn(t, sigma1, sigma2, mu1 = -1, mu2 =1):
    """Plot two gaussian distributions with thresholds t"""
    x = np.linspace(stats.norm.ppf(0.01, loc=mu1, scale=sigma1), stats.norm.ppf(0.99, loc=mu2, scale=sigma2), 100)
    
    plt.figure()
    plt.plot(x, stats.norm.pdf(x, loc=mu1, scale=sigma1), 'black')
    plt.plot(x, stats.norm.pdf(x, loc=mu2, scale=sigma2), 'black')
    for threshold in t:
        plt.plot([threshold, threshold], [0,0.7], '--')

    plt.show()