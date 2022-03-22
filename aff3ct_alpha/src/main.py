#include <aff3ct.hpp>
#using namespace aff3ct;
import sys
import numpy as np

class Source:
        
        def generate(self):
                pass

class Modem:
        
        def modulate(self):
                pass 

        def demodulate(self):
                pass 

class Channel:
        
        def add_noise(self):
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
            "block_length": 18000,
            "storage_pattern": "gray",
            "ditribution_model": "normal",
            "voltage_thresholds": "auto",
            "readout_type": "soft",

            "eb0_min" : 10,
            "eb0_max" : 12,
            "eb0_step": 1
            }

def init_modules(p):
        return {
            "source" : Source(),
            "modem" : Modem(p),
            "channel" : Channel(p)
        }

def init_buffers(p):
        return { 
                "enc_bits" : np.array(),
                "symbols"  : np.array(),
                "noisy_symbols" : np.array(),
        }

def init_utils(m):
        return {}

def main():
        p = init_params()
        m = init_modules(p)
        b = init_buffers(p)
        u = init_utils(m)

#loop over SNRs range
        for ebn0 in np.arange(p["eb0_min"], p["eb0_max"], p["eb0_step"]):
                #compute the current sigma for the channel noise
                esn0  = Tools.ebn0_to_esn0 (ebn0, p["R"]) #TODO
                sigma = Tools.esn0_to_sigma(esn0)         #TODO

                u["noise"].set_values(sigma, ebn0, esn0)

                #update the sigma of the modem and the channel
                m["modem"].set_noise(u["noise"])
                m["channel"].set_noise(u["noise"])

                #display the performance (BER and FER) in real time (in a separate thread)
                #u.terminal->start_temp_report();

                # run the simulation chain
                for i in range(10):
                
                        m["source"].generate(b["ref_bits"])
                        #m.encoder->encode      (b.ref_bits,      b.enc_bits     );
                        m["modem"].modulate(["benc_bits"], b["symbols"])
                        m["channel"].add_noise(b["symbols"], b["noisy_symbols"])
                        m["modem"].demodulate(b["noisy_symbols"], b["LLRs"])
                        #m.decoder->decode_siho (b.LLRs,          b.dec_bits     );
                        #m.monitor->check_errors(b.dec_bits,      b.ref_bits     );
                
                # display the performance (BER and FER) in the terminal
                #u.terminal->final_report();

                #reset the monitor for the next SNR
                #m.monitor->reset();
                #u.terminal->reset();

if __name__ == '__main__':
        main()

