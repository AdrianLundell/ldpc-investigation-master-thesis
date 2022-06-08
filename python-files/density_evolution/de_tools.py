
# %% Imports
import numpy as np
import math
import optimal_thresholds
from scipy.stats import norm


def p_0_symmetrical(probs, min, step, max):
    n_elements = int((max-min)/step)  # Has to be an integer value exactly
    p_0 = np.zeros(n_elements)

    n = probs.size
    for i in range(n):
        llr = math.log(probs[i]/probs[n-1-i])
        idx = math.floor((llr-min)/step)
        p_0[idx] = probs[i]

    return p_0


def mid_point(mu1, mu2, sigma1, sigma2):
    """Calculate the point of intersection for two gaussian distributions"""
    if sigma1 == sigma2:
        return (mu1+mu2)/2

    else:
        a = sigma2**2 - sigma1**2
        b = -(2*sigma2**2*mu1 - 2*sigma1**2*mu2)
        c = sigma2**2*mu1**2 - sigma1**2*mu2**2 - \
            sigma1**2*sigma2**2*np.log(sigma2**2/sigma1**2)
        return (-b + np.sqrt(b**2 - 4*a*c))/(2*a)


def p_yx_and_LLR(t, mu1, mu2, sigma1, sigma2):
    """Calculate LLRs for BI-GAWN channel with thresholds t"""
    p = (np.array([[norm.cdf(t[0], mu1, sigma1), norm.cdf(t[1], mu1, sigma1), norm.cdf(t[2], mu1, sigma1), 1],
                   [norm.cdf(t[0], mu2, sigma2), norm.cdf(t[1], mu2, sigma2), norm.cdf(t[2], mu2, sigma2), 1]]))

    yx_prob = np.array([p[0, 0], p[0, 1] - p[0, 0], p[0, 2] - p[0, 1], p[0, 3] - p[0, 2],
                        p[1, 0], p[1, 1] - p[1, 0], p[1, 2] - p[1, 1], p[1, 3] - p[1, 2]])

    return np.log(yx_prob[0:4]/yx_prob[4:])

# %% Calculating the mutual information between two gaussian distributions


def mutual_info(t):
    p = (np.array([[norm.cdf(t[0], -1, sigma1), norm.cdf(t[1], -1, sigma1), norm.cdf(t[2], -1, sigma1), 1],
                   [norm.cdf(t[0], 1, sigma2), norm.cdf(t[1], 1, sigma2), norm.cdf(t[2], 1, sigma2), 1]]))

    yx_prob = np.array([p[0, 0], p[0, 1] - p[0, 0], p[0, 2] - p[0, 1], p[0, 3] - p[0, 2],
                        p[1, 0], p[1, 1] - p[1, 0], p[1, 2] - p[1, 1], p[1, 3] - p[1, 2]])

    y_prob = np.array([(yx_prob[0]+yx_prob[4])*0.5, (yx_prob[1]+yx_prob[5])*0.5,
                       (yx_prob[2]+yx_prob[6])*0.5, (yx_prob[3]+yx_prob[7])*0.5])

    yx_prob = yx_prob[yx_prob > 0]
    y_prob = y_prob[y_prob > 0]

    y_entropy = -(y_prob @ np.log2(y_prob))
    yx_entropy = -(yx_prob @ np.log2(yx_prob))

    # Note: H(Y|X) = 1/2H(Y|X=0) + 1/2H(Y|X=1)
    # here both terms are concatenated.
    return y_entropy - yx_entropy*0.5


def optimize_threhsolds(mu1, mu2, sigma1, sigma2):
    m = mid_point(mu1, mu2, sigma1, sigma2)
    offsets = np.arange(7.5e-3, 1, 7.5e-3)
    mi_result = np.zeros((len(offsets), len(offsets)))

    for i, positive_offset in enumerate(offsets):
        for j, negative_offset in enumerate(offsets):
            t = [m-negative_offset, m, m+positive_offset]
            mi_result[i, j] = mutual_info(t)

    max_mi = np.max(mi_result)
    max_mi_index = np.unravel_index(np.argmax(mi_result), mi_result.shape)
    pos_offset = offsets[max_mi_index[0]]
    neg_offset = offsets[max_mi_index[1]]

    return max_mi, pos_offset, neg_offset


# %% Calculate many thresholds and output file with columns
#noise_db, sigma_ratio, sigma1, sigma2, max_mi, mid_point, neg_offset, pos_offset, llr1, llr2, llr3, llr4
db_range = np.arange(3, 8, 0.5)
sigma_range = [0.2, 0.3, 0.4, 0.5]

result = np.zeros((len(db_range)*len(sigma_range), 12))
i = 0

for noise_db in db_range:
    for sigma_ratio in sigma_range:
        N0 = 2/10**(noise_db/10)
        sigma1 = np.sqrt(N0*sigma_ratio)
        sigma2 = np.sqrt(N0*(1-sigma_ratio))
        middle = mid_point()

        print(
            f"Calculating thresholds for SNR={noise_db}dB, sigma ratio={sigma_ratio}.")

        max_mi, pos_offset, neg_offset = optimize_threhsolds()
        llr = LLR(np.array([middle-neg_offset, middle, middle+pos_offset]))

        result[i, :] = np.array([noise_db, sigma_ratio, sigma1, sigma2, max_mi,
                                 middle, neg_offset, pos_offset, llr[0], llr[1], llr[2], llr[3]])

        i += 1
