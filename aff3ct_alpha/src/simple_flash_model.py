"""
This script implements the modulation -> channel -> demodulation signal chain for a Flash memory, only modelling neighbouring cells.
Structure is inspired by the AFF3CT library examples for simple integration.
https://aff3ct.readthedocs.io/en/latest/user/library/examples.html
"""
#%%i
from cmath import inf
import numpy as np
import matplotlib.pyplot as plt 
import scipy.stats as sp

class Module():
        """Base class for modules"""
        def __init__(self, p) -> None:
            self.p = p

class Source(Module):
        """Code word generator"""
        def generate(self, ref_bits):
                ref_bits[:] = np.random.randint(2, size=self.p["block_length"])

class Modem(Module):
        """Simulates storage and readout of a FLASH memory cell"""

        def modulate(self, ref_bits, symbols):
                """BPSK modulation"""
                symbols[:] = -ref_bits*2 + 1


        def demodulate(self, noisy_symbols, LLRs):
                """Demodulates a sequence of N/M symbols into a sequence of N*S bits, S being the number of soft bits per read"""
                LLRs[:] = noisy_symbols

class Channel(Module):
        """Simulates all noise sources in the FLASH memory cell"""

        def add_noise(self, symbols, noisy_symbols):
                """BMC derived from average gaussian distributions"""
                T = [-0.4, -0.2, 0, 0.2, 0.4]

                for i in range(len(noisy_symbols)):
                        sigma1 = np.random.rand()
                        sigma2 = np.random.rand()
                        p = np.random.rand()
                        n_bin = 0

                        for t in T:
                                if symbols[i] == -1:
                                        q = sp.norm.cdf(t, -1, sigma1)
                                if symbols[i] == +1:
                                        q = sp.norm.cdf(t, 1, sigma2)
                                if q > p:
                                        break 
                                n_bin += 1

                        noisy_symbols[i] = n_bin


class Monitor(Module):
        def check_errors():
                pass

class Tools:
        
        @staticmethod
        def ebn0_to_esn0():
                pass 

        @staticmethod
        def esn0_to_sigma():
                pass

def init_params():
        return {
            "block_length": 200,            #Length of full encoded bit sequence
            "bit_mapping": "gray",            #Mapping from bit sequences to voltage levels
            "voltage_levels": [0, 2, 4, 6],   #Target levels for each voltage  
            "voltage_distribution" : "AWGN",  #Model of actual voltage distributions
            "threshold_selection": "static_hard",  #Thresholding algorithm
            "threshold_levels" : [-inf, 1, 3, 5, inf],   #Static thresholding levels (inlcude +-inf for lazy implmentation)

            "eb0_min" : 10,
            "eb0_max" : 12,
            "eb0_step": 1
            }

def init_modules(p):
        return {
            "source" : Source(p),
            "modem" : Modem(p),
            "channel" : Channel(p),
            "monitor" : Monitor(p)
        }

def init_buffers(p):
        number_of_symbols = int(p["block_length"]//np.log2(len(p["voltage_levels"])))

        return { 
                "ref_bits" : np.zeros(p["block_length"], dtype=int),
                "symbols"  : np.zeros(p["block_length"]),
                "noisy_symbols" : np.zeros(p["block_length"]),
                "LLRs"  : np.zeros(p["block_length"]),
                "dec_bits"  : np.zeros(p["block_length"]),
                }

def init_utils(m):
        return {}

def main():
        p = init_params()
        m = init_modules(p)
        b = init_buffers(p)
        u = init_utils(m)

#loop over SNRs range
        #for ebn0 in np.arange(p["eb0_min"], p["eb0_max"], p["eb0_step"]):
                #compute the current sigma for the channel noise
                #esn0  = Tools.ebn0_to_esn0 (ebn0, p["R"]) #TODO
                #sigma = Tools.esn0_to_sigma(esn0)         #TODO

                #u["noise"].set_values(sigma, ebn0, esn0)

                #update the sigma of the modem and the channel
                #m["modem"].set_noise(u["noise"])
                #m["channel"].set_noise(u["noise"])

                #display the performance (BER and FER) in real time (in a separate thread)
                #u.terminal->start_temp_report();

                # run the simulation chain
                #for i in range(10):
        
        m["source"].generate(b["ref_bits"])
        #m.encoder->encode      (b.ref_bits,      b.enc_bits     );
        m["modem"].modulate(b["ref_bits"], b["symbols"])
        m["channel"].add_noise(b["symbols"], b["noisy_symbols"])
        m["modem"].demodulate(b["noisy_symbols"], b["LLRs"])

        plt.figure()
        plt.subplot(1,2,1)
        plt.hist(b["noisy_symbols"], 50)
        plt.hist(b["symbols"], 50)                

        #difference = b["ref_bits"] - b["LLRs"]
        #difference = np.reshape(difference, (100,160))
        #plt.subplot(1,2,2)
        #plt.spy(difference)
        
        plt.show()
                #m.decoder->decode_siho (b.LLRs,          b.dec_bits     );
                #m["monitor"].check_errors(b["dec_bits"], b["ref_bits"])
        
                # display the performance (BER and FER) in the terminal
                #u.terminal->final_report();

                #reset the monitor for the next SNR
                #m.monitor->reset();
                #u.terminal->reset();

if __name__ == '__main__':
        main()
