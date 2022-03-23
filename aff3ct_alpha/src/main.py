#include <aff3ct.hpp>
#using namespace aff3ct;
#%%i
from cmath import inf
import numpy as np
import matplotlib.pyplot as plt 

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
                """Modulates a sequence of N bits into a sequence of Q=N/M symbols encdoing M bits per symbol"""
                m = int(np.log2(len(self.p["voltage_levels"])))
                q = int(self.p["block_length"]//m) 
                
                #Interpret bit sequence as a sequence of binary symbols
                #Example: 1001 is rehaped to (10)(01) and converted to (2,1)
                reshaped_bits = np.reshape(ref_bits, (q,m))
                binary_symbols = np.packbits(reshaped_bits, axis=-1, bitorder="little").flatten()

                #Generate mapping from binary symbols to voltage levels
                #Example: (2,1) is in gray code mapped to (3,1) with given voltage levels (v3, v1)
                if self.p["bit_mapping"] == "gray":
                        mapping = {i^(i>>1) : self.p["voltage_levels"][i] for i in range(len(self.p["voltage_levels"]))}
                        
                #Apply mapping
                symbols[:] = np.vectorize(mapping.get)(binary_symbols)
                

        def demodulate(self, noisy_symbols, LLRs):
                """Demodulates a sequence of N/M symbols into a sequence of N*S bits, S being the number of soft bits per read"""
                m = int(np.log2(len(self.p["voltage_levels"])))
                q = int(self.p["block_length"]//m) 
                binary_symbols = np.zeros(q, np.uint8)

                if self.p["threshold_selection"] == "static_hard":
                        if self.p["bit_mapping"] == "gray":
                                mapping = {self.p["threshold_levels"][i+1] : i^(i>>1)for i in range(len(self.p["voltage_levels"]))} 
                        
                        for threshold1, threshold2 in zip(self.p["threshold_levels"][:-1], self.p["threshold_levels"][1:]):
                                binary_symbols[np.logical_and(threshold1 <= noisy_symbols, noisy_symbols < threshold2)] = mapping[threshold2]
                        binary_symbols = np.reshape(binary_symbols, (q,1))

                if self.p["threshold_selection"] == "static_soft":
                        pass 

                if self.p["threshold_selection"] == "dynamic":
                        pass

                reshaped_bits = np.unpackbits(binary_symbols, axis=1, count=m, bitorder="little")
                LLRs[:] = np.reshape(reshaped_bits, self.p["block_length"])

class Channel(Module):
        """Simulates all noise sources in the FLASH memory cell"""

        def add_noise(self, symbols, noisy_symbols):
                """Adds noise to the signal"""
                if self.p["voltage_distribution"] == "AWGN":
                        m = int(np.log2(len(self.p["voltage_levels"])))
                        q = int(self.p["block_length"]//m) 
                        noisy_symbols[:] = symbols + np.random.normal(loc=0, scale=0.0, size=q)


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
            "block_length": 16000,            #Length of full encoded bit sequence
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
                "symbols"  : np.zeros(number_of_symbols),
                "noisy_symbols" : np.zeros(number_of_symbols),
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

        difference = b["ref_bits"] - b["LLRs"]
        difference = np.reshape(difference, (100,160))
        plt.subplot(1,2,2)
        plt.spy(difference)
        
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
