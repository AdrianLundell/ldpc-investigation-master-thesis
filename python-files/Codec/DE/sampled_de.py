import numpy as np 
import matplotlib.pyplot as plt 
import discretize_sigma as d

#%% Test density evolution
def to_cdf(pdf):
    """Returns the discrete pdf of a cdf"""
    return np.cumsum(pdf)

def to_pdf(cdf):
    """Returns the discrete cdf of a pdf"""
    cdf = np.hstack((0, cdf, 1))
    pdf = cdf[1:] - cdf[:-1]
    return pdf[:-1]

def plot_samples(samples, bins, n_samples, cdf = True):
    bins = np.append(bins, np.inf)
    values, bins = np.histogram(samples, bins) 
    values = values/ n_samples
    if cdf: 
        values = to_cdf(values)
    
    plt.plot(bins[:-1], values)

def sampled_density_evolution(pdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, n_samples = 10**4):

    cdf = to_pdf(pdf)

    #Initial messages
    samples = []
    for i in range(n_samples):
        x = np.random.rand()
        sample = f_grid[np.argmax(cdf >= x)]  
        samples.append(sample)
    samples = np.array(samples)

    #gamma
    g0 = -np.log(np.tanh(np.abs(samples[samples<=0])/2))
    g1 = -np.log(np.tanh(np.abs(samples[samples>0])/2))    

    plot_samples(g0, g_grid, n_samples = n_samples)
    plt.show()
    plot_samples(g1, g_grid, n_samples = n_samples)
    plt.show()
    return

    #rho
    sampled_rho0 = []
    for i in range(n_samples//2):
        sample = 0
        for i, coeff in enumerate(rho_coeffs[1:]):
            x = np.random.randint(0, g0.size, i)
            x = np.take(g0, x)
            sample = coeff * np.sum(x)
        sampled_rho0.append(sample)
    sampled_rho0 = np.array(sampled_rho0)

    sampled_rho1 = []
    for i in range(n_samples//2):
        sample = 0
        for i, coeff in enumerate(rho_coeffs[1:]):
            x = np.random.randint(0, g1.size, i)
            x = np.take(g1, x)
            sample = coeff * np.sum(x)
        sampled_rho1.append(sample)
    sampled_rho1 = np.array(sampled_rho1)
    
    #gamma inv
    sampled_inv0 = -np.log((1 + np.exp(-sampled_rho0))/(1 - np.exp(-sampled_rho0)))
    sampled_inv1 = np.log((1 + np.exp(-sampled_rho1))/(1 - np.exp(-sampled_rho1)))
    sampled_inv = np.append(sampled_inv0, sampled_inv1)
    
    #lambda
    sampled_lambda = []
    for i in range(n_samples):
        sample = 0
        for i, coeff in enumerate(lambda_coeffs[1:]):
            x = np.random.randint(0, sampled_inv.size, coeff)
            x = np.take(sampled_inv, x)
            sample = coeff * np.sum(x)
        sampled_lambda.append(sample)
    sampled_lambda = np.array(sampled_lambda)

    #conv 
    sampled_conv = []
    for i in range(n_samples):
        x1 = sampled_inv[np.random.randint(0, n_samples)]
        x2 = samples[np.random.randint(0, n_samples)]
        sample = x1 + x2
        sampled_conv.append(sample)
    sampled_conv = np.array(sampled_conv)



    return sampled_conv

#%%
sigma, p0, bins = d.compute_pdf(0.01, 0.5, 1024, 10)
f_grid, g_grid, pdf = d.create_pdf(p0, bins, 1024)

pl = sampled_density_evolution(pdf, f_grid, g_grid, [0, 0, 0, 0, 0, 1], [0,0,1], 10**6)
#%%
plot_samples(pl, f_grid)
# %%
plt.plot(f_grid, to_cdf(pdf))
# %%
