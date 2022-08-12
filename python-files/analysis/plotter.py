#%%
import rber_to_ber
from config import cfg

if __name__ == '__main__':
    plot_type = cfg.get('plot_type')
    
    if plot_type == "rber_to_ber":
        rber_to_ber.plot()
# %%
