import numpy as np
import matplotlib.pyplot as plt
import plotter_utils as pu
from scipy.stats import norm
import scipy.special as special

from config import cfg


def esn0dB_to_rber(esn0dB):
    esn0 = 10**(esn0dB/10)
    sigma = np.sqrt(1/(2*esn0))
    rber = norm.cdf(0, loc=1, scale=sigma)
    return rber

def plot():
    input_files = cfg.get('input_files') 
    data_dict = pu.read_data(input_files)
    
    for fname in data_dict:
        df = data_dict[fname]
        
        # Convert Es/N0 to RBER
        esn0dB = df["Es/N0"].to_numpy()
        rber = esn0dB_to_rber(esn0dB)
        rberdB = 10*np.log10(rber)

        ber = df["BER"].to_numpy()

        plt.plot(rber, ber, 'o')
    
    plt.legend([fname for fname in data_dict])
    plt.ylabel("BER")
    plt.xlabel('RBER')
    plt.show()
    