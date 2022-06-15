BACKGROUND
Channel model design is an essential part of the ldpc-code optimisation process since the final code is designed to be optimal for the the specific channel. The closer this channel model is to the real case, the more suited the code will be also to the real case.

A NAND flash memory is often modeled as an AWGN-like channel, however this model fails to capture two important characteristics:
    1. Readout is perfromed through applying thresholds, creating a heavily discretized output.
    2. The noise is typically not symmetric 

Because of this a DMC channel with transission probabilities calculated by integrating over the AWGN distributions. The thresholds are calculated by maximising the mutual information of the resulting channel. Furthermore the distributions are not assumed to have the same variance which complicates the calculation of noise levels and thresholds.

Currently only this discretized asymmetric AWGN channel is supported but the ideas should be easily generalizeable to other distributions.

PROGRAM DESCRIPTION
compute_dmc.py is the main program, reading settings from the config.yml. 

Input parameters:
    1. 

Output: A csv-file with the following columns descriving the channel:


SOURCES:
Channel Coding for Nonvolatile Memory Technologies: Theoretical Advances and Practical Considerations:
https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=7920318

Enhanced Precision Through Multiple Reads for LDPC Decoding in Flash Memories: 
https://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=6804933