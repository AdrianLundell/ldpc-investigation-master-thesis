import numpy as np 

def density_evolution(p0_pdf, f_grid, g_grid, n_iter = 50):

    pl = to_cdf(p0_pdf)

    for l in range(n_iter):
        fig,axes = plt.subplots(1,2)
        x1 = gamma(pl, f_grid, g_grid)
        x2 = rho(x1)
        x3 = gamma_inv(x2, f_grid, g_grid)
        x4 = lambd(x3)

        pl = sp.convolve(p0_pdf, x4)
        current_size = pl.size
        pl = pl[:np.argmax(pl)+1]
        pl = np.pad(pl, (0, current_size-pl.size), constant_values = pl.max())

        axes[0].plot(pl)
        axes[1].scatter(l, pl[(2**14*4)//2-1])

        # plt.show()
        # plt.plot(x1[0,:])
        # plt.show()
        # plt.plot(x1[1,:])
        # plt.show()
        # plt.plot(x2[0,:])
        # plt.show()
        # plt.plot(x2[1,:])
        # plt.show()
        # plt.plot(x3)
        # plt.show()

        plt.show()

#%% Test density evolution
def plot_samples(samples, bins, ax, cdf = False):
    bins = np.append(bins, np.inf)
    values, bins = np.histogram(samples, bins) 
    values = values / len(samples)
    if cdf: 
        values = density_to_dist(values)
    ax.plot(bins[:-1], values)

def test(cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, n_samples = 10**4):

    #compare sampled and theoretical cdf
    samples = []
    for i in range(n_samples):
        x = np.random.rand()
        sample = f_grid[np.argmax(cdf >= x)]  
        samples.append(sample)
    samples = np.array(samples)

    fig, ax = plt.subplots(1,1)
    ax.set_title("Compare distributions")
    plot_samples(samples, f_grid, ax, cdf = True)
    plt.plot(f_grid, cdf)
    plt.show()

    #Compare sampled and theoretical gamma(cdf)
    sampled_g0 = -np.log(np.tanh(np.abs(samples[samples<=0])/2))
    sampled_g1 = -np.log(np.tanh(np.abs(samples[samples>0])/2))    
    g = gamma(cdf, f_grid, g_grid)
    
    fig, axes = plt.subplots(1,2, figsize = (15,10))
    fig.suptitle("Compare gamma(distribution)")
    plot_samples(sampled_g0, g_grid, axes[0])
    plot_samples(sampled_g1, g_grid, axes[1])
    axes[0].plot(g_grid, dist_to_density(g[0,:]))
    axes[1].plot(g_grid, dist_to_density(g[1,:]))
    plt.show()

    #Compare sampled and theoretical rho(gamma(cdf)). 
    sampled_rho0 = []
    for i in range(n_samples//2):
        sample = 0
        for i, coeff in enumerate(rho_coeffs[1:]):
            x = np.random.randint(0, sampled_g0.size, i)
            x = np.take(sampled_g0, x)
            sample = coeff * np.sum(x)
        sampled_rho0.append(sample)
    sampled_rho0 = np.array(sampled_rho0)

    sampled_rho1 = []
    for i in range(n_samples//2):
        sample = 0
        for i, coeff in enumerate(rho_coeffs[1:]):
            x = np.random.randint(0, sampled_g1.size, i)
            x = np.take(sampled_g1, x)
            sample = coeff * np.sum(x)
        sampled_rho1.append(sample)
    sampled_rho1 = np.array(sampled_rho1)
    
    fig, axes = plt.subplots(1,2, figsize = (15,10))
    fig.suptitle("Compare rho(gamma(distribution))")
    plot_samples(sampled_rho0, g_grid, axes[0])
    plot_samples(sampled_rho1, g_grid, axes[1])
    plt.show()

    #Compare sampled and theoretical inv gamma
    sampled_inv0 = -np.log((1 + np.exp(-sampled_rho0))/(1 - np.exp(-sampled_rho0)))
    sampled_inv1 = np.log((1 + np.exp(-sampled_rho1))/(1 - np.exp(-sampled_rho1)))
    sampled_inv = np.append(sampled_inv0, sampled_inv1)
    # sampled_pdf, bins = np.histogram(sampled_inv, 30)
    
    # plt.plot(bins[:-1], density_to_dist(sampled_pdf))
    # plt.title("Compare inv_gamma(rho(gamma(distribution)))")
    # plt.show()
    
    #Compare lambda
    sampled_lambda = []
    for i in range(n_samples):
        sample = 0
        for i, coeff in enumerate(lambda_coeffs[1:]):
            x = np.random.randint(0, sampled_inv.size, coeff)
            x = np.take(sampled_inv, x)
            sample = coeff * np.sum(x)
        sampled_lambda.append(sample)
    sampled_lambda = np.array(sampled_lambda)

    #Compare convolution
    sampled_conv = []
    for i in range(n_samples):
        x1 = sampled_inv[np.random.randint(0, n_samples)]
        x2 = samples[np.random.randint(0, n_samples)]
        sample = x1 + x2
        sampled_conv.append(sample)
    sampled_conv = np.array(sampled_conv)

    error = sum(sampled_conv == 0)/n_samples
    print("Error ", error)