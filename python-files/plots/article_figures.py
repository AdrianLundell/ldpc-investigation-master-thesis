from signal import SIGRTMAX
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import sys

sys.path.append("../../compute_dmc")
import dmc_utils as dmc_u

# Valley search algorithm illustration
default_t = 0.15
max_t = 0.5
t_arr = np.linspace(-max_t,max_t,30)

mu = 1
sigma = 0.5
ber = np.zeros(t_arr.size)
for i,t in enumerate(t_arr):
    p0 = norm.cdf(t, -mu, sigma)
    ber_0 = 1-p0
    p1 = norm.cdf(t,mu,sigma)
    ber_1 = p1
    if t < default_t:
        idx_lo = i
    ber[i] = 1/2*(ber_0+ber_1)

t_arr = t_arr-default_t
plt.semilogy(t_arr,ber,marker="s",markersize=2,linewidth=1)
idx = np.argmin(ber)
t_opt = np.full(2,t_arr[idx])
ber_opt = np.array([ber[idx], 0.8*ber[idx]])
plt.semilogy(t_opt, ber_opt, color = '#696969',linestyle="--", linewidth=1)

t_0 = np.full(2, 0)
ber_0 = ber[idx_lo] + (ber[idx_lo+1]-ber[idx_lo])/2
ber_0 = np.array([ber_0, 0.8*ber[idx]])
plt.semilogy(t_0, ber_0, color="k", linestyle="--", linewidth=1)
plt.xlabel("Voltage offset")
plt.ylabel("BER")
plt.show()

# Discrete voltage distributions
n_plots = 3
fig, ax = plt.subplots(n_plots)
titles = ["Single Level Cell (SLC)", "Multi Level Cell (MLC)", "Triple Level Cell (TLC)"]
bit_vals = [['1','0'],['11','10','00','01'],['111','011','001','101','100','000','010','110']]
for n in range(n_plots):
    for i in range(2**(n+1)):
        loc = 6*i
        x = np.linspace(norm.ppf(0.0001, loc=loc), norm.ppf(0.9999,loc=loc), 100)
        ax[n].fill_between(x, norm.pdf(x,loc=loc),0, color='tab:blue')
        if i < 2**(n+1)-1:
            x = np.full(2, (loc+3))
            y = np.array([0,0.5])
            ax[n].plot(x,y,color='k',linestyle='dashed')

        ax[n].text(loc-0.3*n**2, 0.1 , bit_vals[n][i])
        ax[n].set_title(titles[n])
        ax[n].set_xlabel("Voltage level")
        ax[n].set_ylabel("Probability density")
        ax[n].xaxis.set_ticks([])
        ax[n].yaxis.set_ticks([])


fig.tight_layout()
plt.show()

# Mutual information as a function of threshold
t_arr = np.arange(0, 3, 0.001)
mi = np.zeros(t_arr.size)
sigma = 0.5
mu = 1
for i,t in enumerate(t_arr):
    t_list = [-t,0,t]
    mi[i] = dmc_u.mutual_info(t_list,sigma,sigma,-mu,mu)

idx = np.argmax(mi)
t_max = np.full(2, t_arr[idx])
y = np.array([np.min(mi), np.max(mi)])


plt.plot(t_arr,mi)
plt.plot(t_max, y, color='k', linestyle='dashed')
plt.ylabel("Mutual information")
plt.xlabel("Threshold offset, t")
plt.show()


# Hard read
x1 = np.linspace(norm.ppf(0.0001, loc=-mu, scale=sigma),
                 norm.ppf(0.9999, loc=-mu, scale=sigma), 1000)
y1 = norm.pdf(x1, loc=-mu, scale=sigma)
x2 = np.linspace(norm.ppf(0.0001, loc=mu, scale=sigma),
                 norm.ppf(0.9999, loc=mu, scale=sigma), 1000)
y2 = norm.pdf(x1, loc=-mu, scale=sigma)
plt.plot(x1, y1)
plt.plot(x2, y2)
x = np.zeros(2)
y = np.array([0, 1])
plt.fill_between(x1, y1, 0, where = x1>=0, color='r')
plt.fill_between(x2, y2, 0, where = x2<=0, color='r')
plt.plot(x, y, color='k', linestyle='dashed')
plt.ylabel("Probability density")
plt.xlabel("Voltage level")
plt.show()

# Soft read
plt.plot(x1, y1)
plt.plot(x2, y2)
plt.plot(x, y, color='k', linestyle='dashed')

idx = np.argmax(mi)
x_t = np.full(2,t_arr[idx])
plt.plot(x_t, y, color='k', linestyle='dashed')
plt.plot(-x_t, y, color='k', linestyle='dashed')
plt.fill_between(x1, y1, 0, where = x1<=-x_t[0], color='g')
plt.fill_between(x1, y1, 0, where=(x1 >= -x_t[0]) & (x1 <= 0), color='y')
plt.fill_between(x2, y2, 0, where=x2 >= x_t[0], color='c')
plt.fill_between(x2, y2, 0, where = (x2>=0) & (x2<=x_t[0]), color='m')
plt.ylabel("Probability density")
plt.xlabel("Voltage level")
plt.show()



