"""
Analysis of generated noise distributions from Channel_AWGN_asymmetric
"""
from cmath import sqrt
import numpy as np
import matplotlib.pyplot as plt
import math

from config import cfg

# This script assumes that tests in aff3ct have been run with  printData = true in Channel_AWGN_asymmetric_test.cpp
cfg = cfg.get('analyze_randomness')


def plot(data, title, xlabel):
    plt.hist(data, bins=50)
    plt.xlabel(xlabel)
    plt.title(title)
    plt.show()


def compute_average(data):
    return np.mean(data)


def compute_variance(data):
    return np.var(data)


def main():
    sigmas = np.loadtxt(cfg.get('sigmas_file'), dtype=float)
    noise = np.loadtxt(cfg.get('noise_file'), dtype=float)
    params = np.loadtxt(cfg.get('params_file'), dtype=float)
    N = int(params[0])
    sigma_tot = params[1]
    sigma_min = params[2]
    n_sigmas = int(params[3])
    voltage_level = params[4]
    voltage_level_sigma = params[5]

    # Analyze sigmas
    average_sigma = compute_average(sigmas)
    variance_sigma = compute_variance(sigmas)
    sigma_max = sigma_tot + sigma_min*(1-sqrt(sigma_min))
    expected_average_sigma = sigma_min + (sigma_tot-sigma_min)/2
    print("==============================================")
    print("Expected:")
    print("average_sigma = %5.2f" %
          (expected_average_sigma))
    print("Actual:")
    print("average_sigma = %5.2f, variance_sigma : %5.2f" %
          (average_sigma, variance_sigma))
    print
    print("==============================================")
    plot(sigmas, 'Distribution of uniformly generated $\sigma$',
         '$\sigma$')

    # Analyze noise
    average_noise = compute_average(noise)
    variance_noise = compute_variance(noise)
    expected_average = voltage_level
    expected_variance = voltage_level_sigma
    print("==============================================")
    print("Expected:")
    print("average_noise = %5.2f, variance_noise : %5.2f" %
          (expected_average, expected_variance))
    print("Actual:")
    print("average_noise = %5.2f, variance_noise : %5.2f" %
          (average_noise, variance_noise))
    print
    print("==============================================")
    plot(noise, 'Distribution of gaussian noise',
         'voltage')


if __name__ == '__main__':
    main()
