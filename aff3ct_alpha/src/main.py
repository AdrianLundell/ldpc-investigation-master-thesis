#include <aff3ct.hpp>
using namespace aff3ct;
import sys
import numpy as np

def main():
    params = {
            "block_length": 18000,
            "storage_pattern": "gray",
            "ditribution_model": "normal",
            "voltage_thresholds": "auto",
            "readout_type": "soft",
            }
    p = init_params()
    m = init_modules(p)
    b = init_buffers(p)
    u = init_utils(m)

    for ebn0 in np.arange(


{
        params1 p;  init_params1 (p   ); // create and initialize the parameters defined by the user
        modules1 m; init_modules1(p, m); // create and initialize modules
        buffers1 b; init_buffers1(p, b); // create and initialize the buffers required by the modules
        utils1 u;   init_utils1  (m, u); // create and initialize utils

        // display the legend in the terminal
        u.terminal->legend();

        // loop over SNRs range
        for (auto ebn0 = p.ebn0_min; ebn0 < p.ebn0_max; ebn0 += p.ebn0_step)
        {
                // compute the current sigma for the channel noise
                const auto esn0  = tools::ebn0_to_esn0 (ebn0, p.R);
                const auto sigma = tools::esn0_to_sigma(esn0     );

                u.noise->set_values(sigma, ebn0, esn0);

                // update the sigma of the modem and the channel
                m.modem  ->set_noise(*u.noise);
                m.channel->set_noise(*u.noise);

                // display the performance (BER and FER) in real time (in a separate thread)
                u.terminal->start_temp_report();

                // run the simulation chain
                while (!m.monitor->fe_limit_achieved())
                {
                        m.source ->generate    (                 b.ref_bits     );
                        m.encoder->encode      (b.ref_bits,      b.enc_bits     );
                        m.modem  ->modulate    (b.enc_bits,      b.symbols      );
                        m.channel->add_noise   (b.symbols,       b.noisy_symbols);
                        m.modem  ->demodulate  (b.noisy_symbols, b.LLRs         );
                        m.decoder->decode_siho (b.LLRs,          b.dec_bits     );
                        m.monitor->check_errors(b.dec_bits,      b.ref_bits     );
                }

                // display the performance (BER and FER) in the terminal
                u.terminal->final_report();

                // reset the monitor for the next SNR
                m.monitor->reset();
                u.terminal->reset();
        }

        return 0;
}

if __name__ == '__main__':
    main()

