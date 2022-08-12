#%% Imports
import numpy as np
import scipy.stats as stats 
import matplotlib.pyplot as plt 
import dmc_utils 

#%% Preliminries
mu1 = -1
mu2 = 1
sigma = 0.5

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'

#%% Plot disctrete SLC
x = np.linspace(stats.norm.ppf(0.001, loc=mu1, scale=sigma), stats.norm.ppf(0.999, loc=mu2, scale=sigma), 22)
y = np.zeros(22)
y[7] = 0.8
y[15] = 0.8

plt.figure()
plt.step(x, y, 'black')

thresholds = [0]
for threshold in thresholds:
    plt.plot([threshold, threshold], [-0.14,0.8], '--')

plt.text(mu1-0.05, -0.14, "0", size = 20)
plt.text(mu2-0.15, -0.14, "1", size = 20)

plt.axis('off')
plt.ylim((-0.2, 1))
plt.show()

#%% Plot SLC
x = np.linspace(stats.norm.ppf(0.001, loc=mu1, scale=sigma), stats.norm.ppf(0.999, loc=mu2, scale=sigma), 100)

plt.figure()
plt.plot(x, stats.norm.pdf(x, loc=mu1, scale=sigma), 'black')
plt.plot(x, stats.norm.pdf(x, loc=mu2, scale=sigma), 'black')

thresholds = [0]
for threshold in thresholds:
    plt.plot([threshold, threshold], [-0.14,0.8], '--')

plt.text(mu1-0.1, -0.14, "0", size = 20)
plt.text(mu2-0.1, -0.14, "1", size = 20)

plt.axis('off')
plt.ylim((-0.2, 1))
plt.show()


#%% PLot MLC
x1 = np.linspace(stats.norm.ppf(0.001, loc=-3, scale=sigma),-1, 100)
x2 = np.linspace(stats.norm.ppf(0.001, loc=-3, scale=sigma), 1, 100)
x3 = np.linspace(-1, stats.norm.ppf(0.999, loc=3, scale=sigma), 100)
x4 = np.linspace(1, stats.norm.ppf(0.999, loc=3, scale=sigma), 100)

plt.figure()
plt.plot(x1, stats.norm.pdf(x1, loc=-3, scale= sigma), 'black')
plt.plot(x2, stats.norm.pdf(x2, loc=-1, scale=sigma), 'black')
plt.plot(x3, stats.norm.pdf(x3, loc=1, scale=sigma), 'black')
plt.plot(x4, stats.norm.pdf(x4, loc=3, scale=sigma), 'black')

thresholds = [-2, 0, 2]
for threshold in thresholds:
    plt.plot([threshold, threshold], [-0.3,0.8], '--')

plt.text(-3-0.1, -0.14, "0", size = 20)
plt.text(-1-0.1, -0.14, "0", size = 20)
plt.text(1-0.1, -0.14, "1", size = 20)
plt.text(3-0.1, -0.14, "1", size = 20)

plt.text(-3-0.1, -0.26, "0", size = 20)
plt.text(-1-0.1, -0.26, "1", size = 20)
plt.text(1-0.1, -0.26, "1", size = 20)
plt.text(3-0.1, -0.26, "0", size = 20)

plt.axis('off')
plt.ylim((-0.3, 1.2))
plt.show()

#%% Soft information 
x = np.linspace(stats.norm.ppf(0.2, loc=mu1, scale=sigma), stats.norm.ppf(0.8, loc=mu2, scale=sigma), 100)

plt.figure()
plt.plot(x, stats.norm.pdf(x, loc=mu1, scale=sigma), 'black')
plt.plot(x, stats.norm.pdf(x, loc=mu2, scale=sigma), 'black')

thresholds, _ = dmc_utils.optimize_thresholds(sigma, sigma, symmetric=True)
for threshold in thresholds:
    plt.plot([threshold, threshold], [-0.25,0.8], '--')

plt.text(-0.6, -0.14, "0", size = 20)
plt.text(-0.25, -0.14, "0", size = 20)
plt.text(0.15, -0.14, "1", size = 20)
plt.text(0.5, -0.14, "1", size = 20)

plt.text(-0.6, -0.26, "0", size = 20)
plt.text(-0.25, -0.26, "1", size = 20)
plt.text(0.15, -0.26, "1", size = 20)
plt.text(0.5, -0.26, "0", size = 20)

plt.axis('off')
plt.ylim((-0.3, 1))
plt.show()


# %%
