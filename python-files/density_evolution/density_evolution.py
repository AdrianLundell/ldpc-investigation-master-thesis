"""
This file implements the full asymmetric density evolution algorithm for optimisng on the discritized A-BIAWGN channel.

Sources: 
https://arxiv.org/pdf/cs/0509014.pdf
https://arxiv.org/pdf/2001.01249.pdf

"""

from config import cfg
import time
import de_utils as de_u
import numpy as np
import matplotlib.pyplot as plt
import ga_continuous as ga_c
import ga_discrete as ga_d

np.seterr(divide='ignore')
cfg_de = cfg.get('density_evolution')


def main():
    algorithm = cfg_de.get("algorithm")
    if algorithm == "ga_continuous":
        ga_c.ga_continous()

    if algorithm == "ga_discrete":
        ga_d.ga_discrete()


if __name__ == "__main__":
    main()
