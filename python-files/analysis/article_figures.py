import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm



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
        ax[n].set_ylabel("Voltage distribution")
        ax[n].xaxis.set_ticks([])
        ax[n].yaxis.set_ticks([])


fig.tight_layout()
plt.show()
