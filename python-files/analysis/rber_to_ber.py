import numpy as np
import matplotlib.pyplot as plt
import plotter_utils as pu
from scipy.stats import norm
from itertools import cycle
import ntpath

from config import cfg

plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams["figure.figsize"] = (4,5)
plt.rcParams["figure.dpi"] = 300

def esn0dB_to_rber(esn0dB):
    esn0 = 10**(esn0dB/10)
    sigma = np.sqrt(1/(2*esn0))
    rber = norm.cdf(0, loc=1, scale=sigma)
    return rber

def plot():
    input_files = cfg.get('input_files') 
    data_dict = pu.read_data(input_files)
    
    prop_cycle = plt.rcParams['axes.prop_cycle']
    colors = cycle(prop_cycle.by_key()['color'])

    plt.fill_between([1.5e-2, 0.4e-1], [1.e-10, 1.e-10], [1.e0, 1.e0], color="lightgray")
    plt.text(0.3e-2, 0.5e-8, "Error floor")
    plt.text(0.25e-1, 0.5e-4, "Waterfall region")
    
    for fname in data_dict:
        df = data_dict[fname]
        
        # Convert Es/N0 to RBER
        if len(df.columns) == 8:
            esn0dB = df["Es/N0"].to_numpy()
            rber = esn0dB_to_rber(esn0dB)
        elif len(df.columns) == 7:
            rber = df.RBER
            
        ber = df["BER"].to_numpy()
        color = next(colors)
        
        head, tail = ntpath.split(fname)
        label = tail or ntpath.basename(head)
        label = label[:-4]

        plt.plot(rber, ber, c=color, label = label)
        plt.plot(rber, ber, '.', c=color)    

    for limit in cfg.get("limits"):
        color = next(colors)
        
        x = np.array([limit, limit])
        y = np.array([0.1, 1e-9])
        plt.plot(x,y, '--', c=color)

    plt.ylabel("Decoded BER")
    plt.xlabel('RBER')
    plt.grid()
   # plt.legend(loc = 4)
#    plt.xlim([2.5e-3, 1.5e-1])
    plt.xlim([1.5e-3, 1.5e-1])
    plt.ylim([1.e-10, 1.e0])
    plt.xscale("log")
    plt.yscale("log")
    plt.subplots_adjust(bottom=0.15, left = 0.15)
    
    plt.show()
    