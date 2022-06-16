"""
This file implements the full asymmetric density evolution algorithm for optimisng on the discritized A-BIAWGN channel.

Sources: 
https://arxiv.org/pdf/cs/0509014.pdf
https://arxiv.org/pdf/2001.01249.pdf

"""

from config import cfg
import time
import de_utils
import numpy as np
import matplotlib.pyplot as plt
import ga_continuous as gc
import ga_discrete as gd

np.seterr(divide='ignore')
cfg_de = cfg.get('density_evolution')


def main():
    algorithm = cfg_de.get("algorithm")
    print_terminal = cfg_de.get("print_terminal")
    dmc_file_path = cfg_de.get("dmc_file_path")
    n_grid = cfg_de.get("n_grid")

    if algorithm == "ga_continuous":
        gc.ga_continous()

    if algorithm == "ga_discrete":
        gd.ga_discrete()

    def eval(x):
        f_grid, g_grid, pdf = de_utils.init_pdf(
            x, n_grid, file_path=dmc_file_path)
        cdf = de_utils.to_cdf(pdf)
        result = de_utils.symmetric_density_evolution(
            cdf, f_grid, g_grid, rho_coeffs, lambda_coeffs, plot=True)

        return result

    result = de_utils.bisection_search(min, max, eval)


if __name__ == "__main__":
    main()
