## Background

Channel model design is an essential part of the ldpc-code optimisation process since the final code is designed to be optimal for the the specific channel. The closer this channel model is to the real world scenario, the more suited the code should be also in real applications.

A NAND flash memory is often modeled as an AWGN-like channel, however this model fails to capture two important characteristics:
1. Readout is perfromed through applying thresholds, creating a heavily discretized output.
2. The noise is not necessarily symmetric

Because of this a DMC channel is used with transission probabilities calculated by integrating over the AWGN distributions between certain thresholds. The thresholds are calculated by maximising the mutual information of the resulting channel. Furthermore the distributions are not assumed to have the same variance which complicates the calculation of noise levels and thresholds.

Currently only this (skewed) AWGN channel discretized to 1 or 2 bits is supported but the ideas should be easily generalizeable to other distributions. Also, the search algorithm for the mutual information maxima can be improved uppon. 
To use the generated file for noise levels other than the ones used in the output-file interpolation on rber och sigma may be used.


## Program description
compute_dmc.py is the main program, reading settings from the config.yml.

Input parameters:
    n_thresholds: 3                 #3 for soft read, 1 for hard read
    skew: 0.5                       #The distribution of total variance: sigma_tot = sigma1*(1 - skew) + sigma2*skew
    rber_range: [0.001, 0.1]        #Range of RBER
    n_rber:  50                     #Number of rows
    output: "data/4_04_AWGN.csv"    #Name of output-file
    symmetric_thresholds: False     #Optimize threshold offset in both directions from mid-point


Output: 
A csv-file with one row for each noise level and the following columns describing the channel:
    RBER     : Raw bit error-rate of the channel
    SIGMA    : Total variance of the channel
    SIGMA1   : Variance of the first noise distribution
    SIGMA2   : Variance of the second noise distribution
    CAPACITY : The capacity of the channel

## SOURCES:
Channel Coding for Nonvolatile Memory Technologies: Theoretical Advances and Practical Considerations:
https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7920318

Enhanced Precision Through Multiple Reads for LDPC Decoding in Flash Memories: 
https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6804933
