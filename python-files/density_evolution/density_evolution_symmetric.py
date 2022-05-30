
import numpy as np
import matplotlib.pyplot as plt
import math
from pyrsistent import v
from sklearn.feature_selection import f_classif
from numba import types
from numba import jit
from scipy.special import binom
from time import process_time

#from sqlalchemy import false, true
#from sympy import C, factorial2

"""
This script can:
1. Perform density evolution for a fixed lam and rho degree distributions for a symmetric channel
2. Optimize the lam distribution given rho and given a BiAWGN with 3 soft reads
"""

variable_types = [
    ('p_0', types.float64[:]),               # a simple scalar field
    ('p_l', types.float64[:]),
    ('q_l', types.float64[:]),
    ('lam', types.float64[:]),
    ('rho', types.float64[:]),
    ('p_e', types.float64[:]),
    ('max_iter', types.int64),
    ('ext_min', types.float64),
    ('ext_step', types.float64),
    ('ext_size', types.int64),
    ('gamma_min', types.float64),
    ('gamma_step', types.float64),
    ('gamma_size', types.int64)
    # an array field
]


class de_algorithm(object):
    def __init__(self, p_0, ext_min, ext_step, ext_size, gamma_min, gamma_step, gamma_size, lam, rho, max_iter):
        # Input arrays should be numpy arrays
        self.p_0 = p_0
        self.p_l = np.copy(p_0)
        self.q_l = np.zeros(p_0.size)
        self.lam = lam
        self.rho = rho
        self.p_e = np.zeros(max_iter)
        self.max_iter = max_iter
        self.ext_min = ext_min
        self.ext_step = ext_step
        self.ext_size = ext_size
        self.gamma_min = gamma_min
        self.gamma_step = gamma_step
        self.gamma_size = gamma_size

    '''
     def create_dict(self, arr):
        return {
            "min": float(arr[0]),
            "step": float(arr[1]),
            "size": int(arr[2])
        }
   

    @staticmethod
    def binomial(n, m):
        n_ints = np.arange(1, n+1)
        m_ints = n_ints[:m]
        n_m_ints = n_ints[:n-m]
        n_fact = np.prod(n_ints)
        m_fact = np.prod(m_ints)
        n_m_fact = np.prod(n_m_ints)
        return n_fact/(m_fact*n_m_fact)
    '''

    def density_evolution(self):
        for iter in range(self.max_iter):

            # Compute q_l
            self.compute_q_l()

            # Compute p_l
            self.compute_p_l()

            # Error is the integral over probability density of
            # negative messages assuming the all 0 codeword
            self.p_e[iter] = np.sum(
                self.p_l[0:math.ceil(self.ext_size/2)])*self.ext_step
            print(self.p_e[iter])

    def compute_q_l(self):
        # gamma(ext) = (sgn(ext),-lntanh|ext/2|)
        # compute p_gamma(ext) given that we know p_ext=p_l for positive and negative messages
        # Compute p_gamma
        p_gamma_pos, p_gamma_neg, p_gamma_excess = self.ext2gamma()
        prob_zero = self.p_l[math.ceil(self.ext_size/2)]*self.ext_step

        max_dc = self.rho.size
        self.q_l[:] = 0
        for i in range(1, max_dc):  # edges of degree 1 are not valid
            if self.rho[i] != 0:
                convolve_min = self.gamma_min*i
                convolve_step = self.gamma_step
                convolve_size = (self.gamma_size+1)*i-i+1

                f_convolve_pos, f_convolve_neg, prob_result_zero = self.convolve_q_i(
                    p_gamma_pos, p_gamma_neg, p_gamma_excess, i, convolve_min, convolve_step, convolve_size, prob_zero)

                f_n_ext, ofl_pos, ofl_neg = self.gamma2ext(
                    f_convolve_pos, f_convolve_neg, convolve_min, convolve_step, convolve_size, prob_result_zero)

                f_n_ext = self.adjust_for_overflow(
                    f_n_ext, ofl_pos, ofl_neg)

                f_n_ext /= np.sum(f_n_ext)*self.ext_step  # Normalize

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
        ext_pos_min = self.ext_step
        ext_pos_step = self.ext_step
        ext_pos_size = math.floor(self.ext_size/2)
        ext_bin = np.arange(
            ext_pos_min, -self.ext_min+self.ext_step, ext_pos_step)

        bin_lower_bounds = ext_bin - self.ext_step/2
        bin_lower_bounds[0] = ext_bin[0]
        bin_upper_bounds = ext_bin + self.ext_step/2
        bin_upper_bounds[-1] = ext_bin[-1]

        bin_ext_boundaries = np.stack(
            (bin_lower_bounds, bin_upper_bounds), axis=0)
        bin_gamma_boundaries = np.log(np.tanh(bin_ext_boundaries/2))

        # p_gamma_excess = [ofl_pos, ufl_pos, ofl-neg, ufl_neg]
        p_gamma_excess = np.zeros(4)

        # Start with the computation of the positive messages
        p_l_pos = self.p_l[int((self.ext_size-1)/2+1):]
        p_gamma_pos, ofl_pos, ufl_pos = self.map_densities(
            ext_pos_step, ext_pos_size, p_l_pos, bin_gamma_boundaries, self.gamma_min, self.gamma_step, self.gamma_size, "gamma")
        p_gamma_excess[0] = ofl_pos
        p_gamma_excess[1] = ufl_pos

        # Next, the computation of the negative messages
        # Need elements in opposite order because of how we defined our bins
        p_l_neg = np.flip(self.p_l[:int((self.ext_size-1)/2)])
        p_gamma_neg, ofl_neg, ufl_neg = self.map_densities(
            ext_pos_step, ext_pos_size, p_l_neg, bin_gamma_boundaries, self.gamma_min, self.gamma_step, self.gamma_size, "gamma")
        p_gamma_excess[2] = ofl_neg
        p_gamma_excess[3] = ufl_neg

        return p_gamma_pos, p_gamma_neg, p_gamma_excess

    @staticmethod
    @jit
    def map_densities(x_step, x_size, p_x, y_bin_boundaries, y_min, y_step, y_size, y_variable):
        """"
        We know the values of p_x and they are mapped onto the bins represented by y_bin_boundaries.
        Now, we want to map those values to the values that are represented by "y_mapping".
        This is done by thinking thinking of "y_mapping" as a set of bins, where the center points are given by "y_mapping"
        """
        # Centers of the bins
        y_max_idx = y_size-1
        y = np.arange(y_min, y_min+y_step*y_size, y_step)
        p_y = np.zeros(y_size)

        if y_variable == "gamma":
            coeffs = 1/np.sinh(-y)
        elif y_variable == "ext":
            tmp = np.tanh(y/2)
            coeffs = (1/tmp)*((1-tmp**2)*(1/2))

        ofl = 0
        ufl = 0

        # Iterate through each bin and map it onto p_y
        y_min_boundary = y[0]-1/2*y_step
        max_y_boundary = y[-1] + 1/2*y_step

        for i in range(x_size):
            y_idx_range = (y_bin_boundaries[:, i]-y_min)/y_step
            y_idx_range = np.rint(y_idx_range)
            y_lower_idx = int(y_idx_range[0])
            y_upper_idx = int(y_idx_range[1])

            out_of_range = False
            partly_out_of_range = False

            # Entire range lower than the index range of y
            if y_lower_idx < 0 and y_upper_idx < 0:
                ufl += p_x[i]*x_step
                out_of_range = True

            # Entire range above the index range of y
            if y_lower_idx > y_max_idx and y_upper_idx > y_max_idx:
                ofl += p_x[i] * x_step
                out_of_range = True

            # If part of the range is below the index range of y
            if not out_of_range and y_lower_idx < 0:
                y_lower_idx = 0
                y_bin_boundaries[0, i] = y_min_boundary
                partly_out_of_range = True

            # If part of the range is above the index range of y
            if not out_of_range and y_upper_idx > y_max_idx:
                y_upper_idx = y_max_idx
                y_bin_boundaries[1, i] = max_y_boundary
                partly_out_of_range = True

            if not out_of_range:
                # If entire interval maps on a singe index of y
                if y_lower_idx == y_upper_idx:
                    p = p_x[i]*x_step/y_step
                    p_y[y_lower_idx] += p

                # If only maps to two indices of y
                elif y_upper_idx-y_lower_idx == 1:
                    bdy = (y_lower_idx+1/2)*y_step+y_min
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
                    lowbdy = (y_lower_idx+1/2)*y_step + y_min
                    highbdy = (y_upper_idx-1/2)*y_step + y_min
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

                prob_y = np.sum(p_y)*y_step
                prob_x = p_x[i]*x_step
                if prob_y < prob_x:
                    # Underflow (Differs from original script)
                    if y_lower_idx == 0:
                        ufl += prob_x-prob_y
                    else:
                        ofl += prob_x-prob_y

        return p_y, ofl, ufl

    def convolve_q_i(self, p_gamma_pos, p_gamma_neg, p_gamma_excess, i, convolve_min, convolve_step, convolve_size, prob_zero):
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
        size_zeropad = 2**math.ceil(np.log2(size_p*i-i+1))

        # Create zeropadding for message
        zeropad_pos = np.zeros(size_zeropad)
        zeropad_pos[:size_p-1] = p_gamma_pos*convolve_step
        zeropad_pos[size_p] = ofl_pos

        zeropad_neg = np.zeros(size_zeropad)
        zeropad_neg[:size_p-1] = p_gamma_neg*convolve_step
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
        f_result_pos = np.absolute(f_result_pos)/convolve_step
        f_result_neg = np.fft.ifft(F_result_neg)
        f_result_neg = np.absolute(f_result_neg)/convolve_step

        return f_result_pos[:(size_p-1)*i+1], f_result_neg[:(size_p-1)*i+1], prob_result_zero

    def gamma2ext(self, f_gamma_pos, f_gamma_neg, convolve_min, convolve_step, convolve_size, prob_result_zero):
        """
        We split p_l in positive and negative parts to preserve information about the sign
        of messages.
        """
        # We do not include the zero-point in our interval since we are dealing with that separately
        ext_pos_min = self.ext_step
        ext_pos_step = self.ext_step
        ext_pos_size = math.floor(self.ext_size/2)

        bin_gamma = np.arange(convolve_min, convolve_min +
                              convolve_size*convolve_step, convolve_step)
        bin_gamma_lower_bounds = bin_gamma - convolve_step/2
        bin_gamma_lower_bounds[0] = bin_gamma[0]
        bin_gamma_upper_bounds = bin_gamma + convolve_step/2
        bin_gamma_upper_bounds[-1] = bin_gamma[-1]
        bin_gamma_boundaries = np.stack(
            (bin_gamma_lower_bounds, bin_gamma_upper_bounds), axis=0)
        ext_bin_boundaries = 2*np.arctanh(np.exp(bin_gamma_boundaries))

        f_ext_pos, ofl_pos, ufl_pos = self.map_densities(
            convolve_step, convolve_size, f_gamma_pos, ext_bin_boundaries, ext_pos_min, ext_pos_step, ext_pos_size, "ext")
        f_ext_neg, ofl_neg, ufl_neg = self.map_densities(
            convolve_step, convolve_size, f_gamma_neg, ext_bin_boundaries, ext_pos_min, ext_pos_step, ext_pos_size, "ext")

        # Flip the negative part back
        f_ext_neg = np.flip(f_ext_neg)

        # Compute the density at zero
        f_ext_0 = (prob_result_zero + ufl_pos + ufl_neg)/self.ext_step
        f_ext = np.concatenate((f_ext_neg, [f_ext_0], f_ext_pos))

        return f_ext, ofl_pos, ofl_neg

    def adjust_for_overflow(self, f_ext, ofl_pos, ofl_neg):

        ext_step = self.ext_step
        f_ext[0] += ofl_neg/ext_step
        f_ext[-1] += ofl_pos/ext_step

        return f_ext

    def compute_p_l(self):

        max_dv = self.lam.size
        self.p_l[:] = 0
        min_ext = self.ext_min
        ext_step = self.ext_step
        ext_elements = self.ext_size
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
