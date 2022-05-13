# %% Help functions
import numpy as np
import matplotlib.pyplot as plt
import math
from pyrsistent import v
from scipy.special import binom
from sklearn.feature_selection import f_classif

from sqlalchemy import false, true
from sympy import C, factorial2

"""
This script can:
1. Perform density evolution for a fixed lam and rho degree distributions for a symmetric channel
2. Optimize the lam distribution given rho and given a BiAWGN with 3 soft reads
"""

# %%


class de_algorithm:
    def __init__(self, p_0, ext, gamma, lam, rho, max_iter):
        # Input arrays should be numpy arrays
        self.p_0 = p_0
        self.p_l = np.copy(p_0)
        self.q_l = np.zeros(p_0.size)
        self.ext = self.create_dict(ext)
        self.gamma = self.create_dict(gamma)
        self.lam = lam
        self.rho = rho
        self.max_iter = max_iter
        self.p_e = np.zeros(max_iter)

    def create_dict(self, arr):
        return {
            "min": float(arr[0]),
            "step": float(arr[1]),
            "size": int(arr[2])
        }

    def density_evolution(self):

        for iter in range(self.max_iter):

            # Compute q_l
            self.compute_q_l()

            # Compute p_l
            self.compute_p_l()

            # Error is the integral over probability density of
            # negative messages assuming the all 0 codeword
            self.p_e[iter] = np.sum(
                self.p_l[0:math.ceil(self.ext['size']/2)])*self.ext['step']

    def compute_q_l(self):
        # gamma(ext) = (sgn(ext),-lntanh|ext/2|)
        # compute p_gamma(ext) given that we know p_ext=p_l for positive and negative messages
        # Compute p_gamma
        p_gamma_pos, p_gamma_neg, p_gamma_excess = self.ext2gamma()
        prob_zero = self.p_l[math.ceil(self.ext['size']/2)]*self.ext['step']

        max_dc = self.rho.size
        self.q_l[:] = 0
        for i in range(1, max_dc):  # edges of degree 1 are not valid
            if self.rho[i] != 0:
                n_log = self.create_dict([self.gamma['min']*i, self.gamma['step'],
                                          (self.gamma['size']+1)*i-i+1])

                f_n_log_pos, f_n_log_neg, prob_result_zero = self.convolve_q_i(
                    p_gamma_pos, p_gamma_neg, p_gamma_excess, i, n_log, prob_zero)

                f_n_ext, ofl_pos, ofl_neg = self.gamma2ext(
                    f_n_log_pos, f_n_log_neg, n_log, prob_result_zero)

                f_n_ext = self.adjust_for_overflow(
                    f_n_ext, ofl_pos, ofl_neg)

                f_n_ext /= np.sum(f_n_ext)*self.ext['step']  # Normalize

                self.q_l += self.rho[i]*f_n_ext
        # gamma(ext) = (sgn(ext),-lntanh|ext/2|)
        # compute p_gamma(ext) given that we know p_ext=p_l for positive and negative messages

        # Have to treat the value of p_l(0) seperately

    def ext2gamma(self):
        """
        We split p_l in positive and negative parts to preserve information about the sign
        of messages.
        """
        # We do not include the zero-point in our interval since we are dealing with that separately
        ext_bin_generator = self.create_dict([
            self.ext['step'], self.ext['step'], math.floor(self.ext['size']/2)])
        ext_bin = np.arange(
            self.ext['step'], -self.ext['min']+self.ext['step'], self.ext['step'])

        bin_lower_bounds = ext_bin - self.ext['step']/2
        bin_lower_bounds[0] = ext_bin[0]
        bin_upper_bounds = ext_bin + self.ext['step']/2
        bin_upper_bounds[-1] = ext_bin[-1]

        bin_ext_boundaries = np.stack(
            (bin_lower_bounds, bin_upper_bounds), axis=0)
        bin_gamma_boundaries = np.log(np.tanh(bin_ext_boundaries/2))

        # p_gamma_excess = [ofl_pos, ufl_pos, ofl-neg, ufl_neg]
        p_gamma_excess = np.zeros(4)

        # Start with the computation of the positive messages
        p_l_pos = self.p_l[int((self.ext['size']-1)/2+1):]
        p_gamma_pos, ofl_pos, ufl_pos = self.map_densities(
            ext_bin_generator, p_l_pos, bin_gamma_boundaries, self.gamma, "gamma")
        p_gamma_excess[0] = ofl_pos
        p_gamma_excess[1] = ufl_pos

        # Next, the computation of the negative messages
        # Need elements in opposite order because of how we defined our bins
        p_l_neg = np.flip(self.p_l[:int((self.ext['size']-1)/2)])
        p_gamma_neg, ofl_neg, ufl_neg = self.map_densities(
            ext_bin_generator, p_l_neg, bin_gamma_boundaries, self.gamma, "gamma")
        p_gamma_excess[2] = ofl_neg
        p_gamma_excess[3] = ufl_neg

        return p_gamma_pos, p_gamma_neg, p_gamma_excess

    def map_densities(self, x_generator, p_x, y_bin_boundaries, y_mapping, y_variable):
        """"
        We know the values of p_x and they are mapped onto the bins represented by y_bin_boundaries.
        Now, we want to map those values to the values that are represented by "y_mapping".
        This is done by thinking thinking of "y_mapping" as a set of bins, where the center points are given by "y_mapping"
        """
        min_y = y_mapping['min']
        y_step = y_mapping['step']
        max_y_idx = y_mapping['size']-1
        # Centers of the bins
        y = np.arange(min_y, min_y+y_step*y_mapping['size'], y_step)
        p_y = np.zeros(max_y_idx+1)

        if y_variable == "gamma":
            coeffs = 1/np.sinh(-y)
        elif y_variable == "ext":
            tmp = np.tanh(y/2)
            coeffs = (1/tmp)*((1-tmp**2)*(1/2))

        ofl = 0
        ufl = 0
        ext_step = x_generator['step']

        # Iterate through each bin and map it onto p_y
        min_y_boundary = y[0]-1/2*y_step
        max_y_boundary = y[-1] + 1/2*y_step
        for i in range(x_generator['size']):
            y_idx_range = np.round(
                (y_bin_boundaries[:, i]-y_mapping['min'])/y_mapping['step']).astype(int)
            y_lower_idx = y_idx_range[0]
            y_upper_idx = y_idx_range[1]

            out_of_range = False
            partly_out_of_range = False

            # Entire range lower than the index range of y
            if y_lower_idx < 0 and y_upper_idx < 0:
                ufl += p_x[i]*ext_step
                out_of_range = true

            # Entire range above the index range of y
            if y_lower_idx > max_y_idx and y_upper_idx > max_y_idx:
                ofl += p_x[i] * ext_step
                out_of_range = True

            # If part of the range is below the index range of y
            if not out_of_range and y_lower_idx < 0:
                y_lower_idx = 0
                y_bin_boundaries[0, i] = min_y_boundary
                partly_out_of_range = True

            # If part of the range is above the index range of y
            if not out_of_range and y_upper_idx > max_y_idx:
                y_upper_idx = max_y_idx
                y_bin_boundaries[1, i] = max_y_boundary
                partly_out_of_range = True

            if not out_of_range:
                # If entire interval maps on a singe index of y
                if y_lower_idx == y_upper_idx:
                    p = p_x[i]*ext_step/y_step
                    p_y[y_lower_idx] += p

                # If only maps to two indices of y
                elif y_upper_idx-y_lower_idx == 1:
                    bdy = (y_lower_idx+1/2)*y_step+min_y
                    lowfrac = abs(y_bin_boundaries[0, i] - bdy)/y_step
                    highfrac = abs(y_bin_boundaries[1, i]-bdy)/y_step
                    p = np.zeros(2)
                    p[0] += coeffs[y_lower_idx] * \
                        lowfrac*p_x[i]
                    p[1] += coeffs[y_upper_idx] * \
                        highfrac*p_x[i]
                    p_y[y_lower_idx: y_upper_idx+1] += p

                # If maps onto multiple indices of y
                else:
                    lowbdy = (y_lower_idx+1/2)*y_step + min_y
                    highbdy = (y_upper_idx-1/2)*y_step + min_y
                    lowfrac = abs(y_bin_boundaries[0, i]-lowbdy)/y_step
                    highfrac = abs(y_bin_boundaries[1, i]-highbdy)/y_step

                    p = np.zeros(y_upper_idx-y_lower_idx+1)
                    p[0] = coeffs[y_lower_idx] * \
                        lowfrac*p_x[i]
                    p[-1] = coeffs[y_upper_idx] * \
                        highfrac*p_x[i]
                    p[1:-1] = coeffs[y_lower_idx+1:y_upper_idx] * \
                        p_x[i]  # Ser skumt ut
                    p_y[y_lower_idx:y_upper_idx+1] += p

            if partly_out_of_range:

                prob_y = np.sum(p)*y_step
                prob_ext = p_x[i]*ext_step
                if prob_y < prob_ext:
                    # Underflow (Differs from original script)
                    if y_lower_idx == 1:
                        ufl += prob_ext-prob_y
                    if y_upper_idx == max_y_idx:
                        ofl += prob_ext-prob_y

        return p_y, ofl, ufl

    def map_densities_fast(self, x_generator, p_x, y_bin_boundaries, y_mapping, y_variable):
        """"
        We know the values of p_x and they are mapped onto the bins represented by y_bin_boundaries.
        Now, we want to map those values to the values that are represented by "y_mapping".
        This is done by thinking thinking of "y_mapping" as a set of bins, where the center points are given by "y_mapping"
        """
        min_y = y_mapping['min']
        y_step = y_mapping['step']
        max_y_idx = y_mapping['size']-1
        # Centers of the bins
        y = np.arange(min_y, min_y+y_step*y_mapping['size'], y_step)
        p_y = np.zeros(max_y_idx+1)

        if y_variable == "gamma":
            coeffs = 1/np.sinh(-y)
        elif y_variable == "ext":
            tmp = np.tanh(y/2)
            coeffs = (1/tmp)*((1-tmp**2)*(1/2))

        ofl = 0
        ufl = 0
        ext_step = x_generator['step']

        # Iterate through each bin and map it onto p_y
        min_y_boundary = y[0]-1/2*y_step
        max_y_boundary = y[-1] + 1/2*y_step
        for i in range(x_generator['size']):
            y_idx_range = np.round(
                (y_bin_boundaries[:, i]-y_mapping['min'])/y_mapping['step']).astype(int)
            y_lower_idx = y_idx_range[0]
            y_upper_idx = y_idx_range[1]

            out_of_range = False
            partly_out_of_range = False

            # Entire range lower than the index range of y
            if y_lower_idx < 0 and y_upper_idx < 0:
                ufl += p_x[i]*ext_step
                out_of_range = true

            # Entire range above the index range of y
            if y_lower_idx > max_y_idx and y_upper_idx > max_y_idx:
                ofl += p_x[i] * ext_step
                out_of_range = True

            # If part of the range is below the index range of y
            if not out_of_range and y_lower_idx < 0:
                y_lower_idx = 0
                y_bin_boundaries[0, i] = min_y_boundary
                partly_out_of_range = True

            # If part of the range is above the index range of y
            if not out_of_range and y_upper_idx > max_y_idx:
                y_upper_idx = max_y_idx
                y_bin_boundaries[1, i] = max_y_boundary
                partly_out_of_range = True

            if not out_of_range:
                # If entire interval maps on a singe index of y
                if y_lower_idx == y_upper_idx:
                    p = p_x[i]*ext_step/y_step
                    p_y[y_lower_idx] += p

                # If only maps to two indices of y
                elif y_upper_idx-y_lower_idx == 1:
                    bdy = (y_lower_idx+1/2)*y_step+min_y
                    lowfrac = abs(y_bin_boundaries[0, i] - bdy)/y_step
                    highfrac = abs(y_bin_boundaries[1, i]-bdy)/y_step
                    p = np.zeros(2)
                    p[0] += coeffs[y_lower_idx] * \
                        lowfrac*p_x[i]
                    p[1] += coeffs[y_upper_idx] * \
                        highfrac*p_x[i]
                    p_y[y_lower_idx: y_upper_idx+1] += p

                # If maps onto multiple indices of y
                else:
                    lowbdy = (y_lower_idx+1/2)*y_step + min_y
                    highbdy = (y_upper_idx-1/2)*y_step + min_y
                    lowfrac = abs(y_bin_boundaries[0, i]-lowbdy)/y_step
                    highfrac = abs(y_bin_boundaries[1, i]-highbdy)/y_step

                    p = np.zeros(y_upper_idx-y_lower_idx+1)
                    p[0] = coeffs[y_lower_idx] * \
                        lowfrac*p_x[i]
                    p[-1] = coeffs[y_upper_idx] * \
                        highfrac*p_x[i]
                    p[1:-1] = coeffs[y_lower_idx+1:y_upper_idx] * \
                        p_x[i]  # Ser skumt ut
                    p_y[y_lower_idx:y_upper_idx+1] += p

            if partly_out_of_range:

                prob_y = np.sum(p)*y_step
                prob_ext = p_x[i]*ext_step
                if prob_y < prob_ext:
                    # Underflow (Differs from original script)
                    if y_lower_idx == 1:
                        ufl += prob_ext-prob_y
                    if y_upper_idx == max_y_idx:
                        ofl += prob_ext-prob_y

        return p_y, ofl, ufl

    def convolve_q_i(self, p_gamma_pos, p_gamma_neg, p_gamma_excess, i, n_log, prob_zero):
        """
        In order to do the convolution, we compute the FFTs of gamma probabilities.
        Thus, the convolution results in a multiplication of the FFTs.

        ofl_pos is the probability of a gamma message being zero. This is directly included
        in the sum. Gamma can only be zero if the sign is positive
        """

        ofl_pos = p_gamma_excess[0]
        ufl_pos = p_gamma_excess[1]
        ofl_neg = p_gamma_excess[2]
        ufl_neg = p_gamma_excess[3]

        # Size of messages, including the zero message corresponding to ofl
        size_p = p_gamma_pos.size + 1
        size_zeropad = 2**math.ceil(math.log2(size_p*i-i+1))

        # Create zeropadding for message
        zeropad_pos = np.zeros(size_zeropad)
        zeropad_pos[:size_p-1] = p_gamma_pos*n_log['step']
        zeropad_pos[size_p] = ofl_pos

        zeropad_neg = np.zeros(size_zeropad)
        zeropad_neg[:size_p-1] = p_gamma_neg*n_log['step']
        zeropad_neg[size_p] = ofl_neg

        prob_pos_fin = np.sum(zeropad_pos)
        prob_pos = prob_pos_fin + ufl_pos

        prob_neg_fin = np.sum(zeropad_neg)
        prob_neg = prob_neg_fin + ufl_neg

        F_pos_fin = np.fft.fft(zeropad_pos/prob_pos_fin)
        F_neg_fin = np.fft.fft(zeropad_neg/prob_neg_fin)

        F_result_pos = np.zeros(size_zeropad, dtype=np.complex)
        F_result_neg = np.zeros(size_zeropad, dtype=np.complex)
        prob_result_zero = 0

        for c in range(i+1):
            prob_c_zero = (1-(prob_pos/(prob_pos+prob_zero)))**(i-c) * \
                binom(i, c)*prob_neg**c*(prob_pos+prob_zero)**(i-c)
            prob_result_zero += prob_c_zero

            factor1 = (1-ufl_neg)**c*(1-ufl_pos)**(i-c)
            factor2 = binom(i, c)*prob_neg**c*prob_pos**(i-c)
            F_c = F_pos_fin**(i-c)*F_neg_fin**c
            F_c *= factor1*factor2

            if c % 2 == 0:
                F_result_pos += F_c
            else:
                F_result_neg += F_c

            for d in range(i-c+1):
                w = (1-ufl_pos)**c*(1-ufl_neg)**d
                v = (1-w)*math.factorial(i)/(math.factorial(c)
                                             * math.factorial(d)*math.factorial(i-c-d))
                v = v*prob_neg**c*prob_pos**d*prob_zero**(i-c-d)
                prob_result_zero += v

        f_result_pos = np.fft.ifft(F_result_pos)
        f_result_pos = np.absolute(f_result_pos)/n_log['step']
        f_result_neg = np.fft.ifft(F_result_neg)
        f_result_neg = np.absolute(f_result_neg)/n_log['step']

        return f_result_pos[:(size_p-1)*i+1], f_result_neg[:(size_p-1)*i+1], prob_result_zero

    def gamma2ext(self, f_n_log_pos, f_n_log_neg, n_log, prob_result_zero):
        """
        We split p_l in positive and negative parts to preserve information about the sign
        of messages.
        """
        # We do not include the zero-point in our interval since we are dealing with that separately
        ext_mapping = self.create_dict([
            self.ext['step'], self.ext['step'], math.ceil((self.ext['size']-1)/2)])

        bin_gamma = np.arange(n_log['min'], n_log['min'] +
                              n_log['size']*n_log['step'], n_log['step'])
        bin_gamma_lower_bounds = bin_gamma - n_log['step']/2
        bin_gamma_lower_bounds[0] = bin_gamma[0]
        bin_gamma_upper_bounds = bin_gamma + n_log['step']/2
        bin_gamma_upper_bounds[-1] = bin_gamma[-1]
        bin_gamma_boundaries = np.stack(
            (bin_gamma_lower_bounds, bin_gamma_upper_bounds), axis=0)
        ext_bin_boundaries = 2*np.arctanh(np.exp(bin_gamma_boundaries))

        f_ext_pos, ofl_pos, ufl_pos = self.map_densities(
            n_log, f_n_log_pos, ext_bin_boundaries, ext_mapping, "ext")
        f_ext_neg, ofl_neg, ufl_neg = self.map_densities(
            n_log, f_n_log_neg, ext_bin_boundaries, ext_mapping, "ext")

        # Flip the negative part back
        f_ext_neg = np.flip(f_ext_neg)

        # Compute the density at zero
        f_ext_0 = (prob_result_zero + ufl_pos + ufl_neg)/self.ext['step']
        f_ext = np.concatenate((f_ext_neg, [f_ext_0], f_ext_pos))

        return f_ext, ofl_pos, ofl_neg

    def adjust_for_overflow(self, f_ext, ofl_pos, ofl_neg):

        ext_step = self.ext['step']
        f_ext[0] += ofl_neg/ext_step
        f_ext[-1] += ofl_pos/ext_step

        return f_ext

    def compute_p_l(self):

        max_dv = self.lam.size
        self.p_l[:] = 0
        min_ext = self.ext['min']
        ext_step = self.ext['step']
        ext_elements = self.ext['size']
        size_p_0 = self.p_0.size
        size_q_l = self.q_l.size
        for j in range(1, max_dv):  # Edges of degree 1 are not valid
            if self.lam[j] != 0:
                size_zeropad = size_p_0+size_q_l*j-j
                zeropad_p_0 = np.zeros(size_zeropad)
                zeropad_p_0[:size_p_0] = self.p_0*ext_step

                zeropad_q = np.zeros(size_zeropad)
                zeropad_q[:size_q_l] = self.q_l*ext_step

                F_p_0 = np.fft.fft(zeropad_p_0)
                F_q = np.fft.fft(zeropad_q)**j

                F_p = F_p_0*F_q

                zeropad_p = np.fft.ifft(F_p)

                minx = min_ext + min_ext*j

                ext_min_idx = round((min_ext-minx)/ext_step)
                ext_max_idx = ext_min_idx+ext_elements-1

                ufl = abs(np.sum(zeropad_p[:ext_min_idx]))
                ofl = abs(np.sum(zeropad_p[ext_max_idx+1:]))

                f_ext = zeropad_p[ext_min_idx:ext_max_idx+1]
                f_ext[0] = ufl
                f_ext[-1] = ofl

                f_ext = np.absolute(f_ext)/ext_step

                self.p_l += self.lam[j]*f_ext


# %%
# Degree distributions
lam = np.array([0, 0.2895, 0.3158, 0, 0, 0.3947])
rho = np.array([0, 0, 0, 0, 0, 0.9032, 0.0968])

# Sample channel (probability density of messages)
# Represents ext = np.arange(-30,30.01, 0.01)

# Extrinsic messages generator [min,step,size]
ext = [-30, 0.01, 6001]
p_0 = np.zeros(6001)
p_0[3239] = 91.606
p_0[2761] = 8.394
# Notice that the integral over the channel becomes 1

# Extrinsic messages are mapped onto gamma [min,step,size]
gamma = [-10, 0.0002, 50000]

# Density evolution
max_iter = 1

alg = de_algorithm(p_0, ext, gamma, lam, rho, max_iter)
alg.density_evolution()

print(alg.p_e)
# %%
