#include <iostream>
#include <memory>
#include <vector>
#include <string>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>


using namespace aff3ct;

struct params 
{
	int K = 16384;		         // number of information bits
	int N = 18592;	 		     // codeword size
	int fe = 100;			 // number of frame errors
	int seed = 0;			 // PRNG seed for the AWGN channel
	float ebn0_min = .5f;	 // minimum SNR value
	float ebn0_max = 2.00f;  // maximum SNR value
	float ebn0_step = .5f;   // SNR step
	float R;				 // code rate (R=K/N)
	float min_sigma = 0.05;  // Set minimum value for individual sigmas
	
    std::vector<float> voltage_levels{-1, 1};
    std::string reader_fpath = "reader_static_hard.txt";
    int page_type = tools::Flash_reader<float, float>::lower;
    int read_type = tools::Flash_reader<float, float>::hard;
    int cell_type = tools::Flash_cell::SLC;

	// std::vector<uint32_t> info_bits = {0,1,2,3,4,5,6,7,8,9,10,11};

	std::vector<uint32_t> info_bits = std::vector<uint32_t>(K);
};
void init_params(params &p);

struct modules
{
	std::unique_ptr<module::Source_random<>> source;
	std::unique_ptr<module::Encoder_LDPC_from_QC<>> encoder;
	std::unique_ptr<module::Modem_flash_page<>> modem;
	std::unique_ptr<module::Channel_AWGN_asymmetric<>> channel;
	std::unique_ptr<module::Decoder_LDPC_BP_flooding_SPA<>> decoder;
	std::unique_ptr<module::Monitor_BFER<>> monitor;
};
void init_modules(const params &p, modules &m);

struct buffers
{
	std::vector<int> ref_bits;
	std::vector<int> enc_bits;
	std::vector<unsigned> voltage_level_indexes;
	std::vector<float> noisy_voltage_levels;
	std::vector<float> noisy_bits;
	std::vector<int> dec_bits;
};
void init_buffers(const params &p, buffers &b);

struct utils
{
	std::unique_ptr<tools::Sigma_asymmetric<>> noise;	     // a sigma noise type
	std::vector<std::unique_ptr<tools::Reporter>> reporters; // list of reporters dispayed in the terminal
	std::unique_ptr<tools::Terminal_std> terminal;			 // manage the output text in the terminal
};
void init_utils(const modules &m, utils &u);

int main(int argc, char **argv)
{
	// get the AFF3CT version
	const std::string v = "v" + std::to_string(tools::version_major()) + "." +
						  std::to_string(tools::version_minor()) + "." +
						  std::to_string(tools::version_release());

	std::cout << "#----------------------------------------------------------" << std::endl;
	std::cout << "# This is a basic program using the AFF3CT library (" << v << ")" << std::endl;
	std::cout << "# Feel free to improve it as you want to fit your needs." << std::endl;
	std::cout << "#----------------------------------------------------------" << std::endl;
	std::cout << "#" << std::endl;

	params p;
	init_params(p); // create and initialize the parameters defined by the user
	modules m;
	init_modules(p, m); // create and initialize the modules
	buffers b;
	init_buffers(p, b); // create and initialize the buffers required by the modules
	utils u;
	init_utils(m, u); // create and initialize the utils

	// display the legend in the terminal
	u.terminal->legend();

	// loop over the various SNRs
	for (auto ebn0 = p.ebn0_min; ebn0 < p.ebn0_max; ebn0 += p.ebn0_step)
	{
		// compute the current sigma for the channel noise
		const auto esn0 = tools::ebn0_to_esn0(ebn0, p.R);
		const auto tot_sigma = tools::esn0_to_sigma(esn0);

		//u.noise->set_noise(sigma, ebn0, esn0);

		u.noise->set_sigmas(tot_sigma, p.voltage_levels.size(), p.min_sigma ,ebn0, esn0);
		// update the sigma of the modem and the channel
		m.channel->set_noise(*u.noise);

		// display the performance (BER and FER) in real time (in a separate thread)
		u.terminal->start_temp_report();

		// run the simulation chain
		while (!m.monitor->fe_limit_achieved() && !u.terminal->is_interrupt())
		{
			m.source->generate(b.ref_bits);
            m.encoder->encode(b.ref_bits, b.enc_bits);
            m.modem->modulate(b.enc_bits, b.voltage_level_indexes);          
     
	        m.channel->generate_sigmas();
	        m.modem->set_noise(*m.channel);
			
			m.channel->add_noise(b.voltage_level_indexes.data(), b.noisy_voltage_levels.data());
			m.modem->demodulate(b.noisy_voltage_levels, b.noisy_bits);
			m.decoder->decode_siho(b.noisy_bits, b.dec_bits);
			m.monitor->check_errors(b.dec_bits, b.ref_bits);
		}

		// display the performance (BER and FER) in the terminal
		u.terminal->final_report();

		// reset the monitor for the next SNR
		m.monitor->reset();
		u.terminal->reset();

		// if user pressed Ctrl+c twice, exit the SNRs loop
		if (u.terminal->is_over())
			break;
	}

	std::cout << "# End of the simulation" << std::endl;

	return 0;
}

void init_params(params &p)
{
	p.R = (float)p.K / (float)p.N;

	for (auto i = 0; i < p.K; i++)
		p.info_bits[i] = i;

	std::cout << "# * Simulation parameters: " << std::endl;
	std::cout << "#    ** Frame errors   = " << p.fe << std::endl;
	std::cout << "#    ** Noise seed     = " << p.seed << std::endl;
	std::cout << "#    ** Info. bits (K) = " << p.K << std::endl;
	std::cout << "#    ** Frame size (N) = " << p.N << std::endl;
	std::cout << "#    ** Code rate  (R) = " << p.R << std::endl;
	std::cout << "#    ** SNR min   (dB) = " << p.ebn0_min << std::endl;
	std::cout << "#    ** SNR max   (dB) = " << p.ebn0_max << std::endl;
	std::cout << "#    ** SNR step  (dB) = " << p.ebn0_step << std::endl;
	std::cout << "#" << std::endl;
}


void init_modules(const params &p, modules &m)
{
	m.source = std::unique_ptr<module::Source_random<>>(new module::Source_random<>(p.K));
	
	const tools::Sparse_matrix H = tools::LDPC_matrix_handler::read("qc_test.qc");
	m.encoder = std::unique_ptr<module::Encoder_LDPC_from_QC<>>(new module::Encoder_LDPC_from_QC<>(p.K, p.N, H));
 
    tools::Flash_cell cell(p.cell_type);
    tools::Flash_reader<float, float> reader(p.page_type, p.read_type, p.reader_fpath);
	m.modem = std::unique_ptr<module::Modem_flash_page<>>(new module::Modem_flash_page<>(p.N, cell, reader));

    tools::Sigma_asymmetric<float> sigmas;
    sigmas.set_sigmas(1, 2, 0.1 , 0, 0);

	m.channel = std::unique_ptr<module::Channel_AWGN_asymmetric<>>(new module::Channel_AWGN_asymmetric<float>(p.N, p.voltage_levels, sigmas));
	
    m.decoder = std::unique_ptr<module::Decoder_LDPC_BP_flooding_SPA<>>(new module::Decoder_LDPC_BP_flooding_SPA<>(
	p.K, 
	p.N, 
	10, 
	H,
	p.info_bits));
	m.monitor = std::unique_ptr<module::Monitor_BFER<>>(new module::Monitor_BFER<>(p.K, p.fe));
};

void init_buffers(const params &p, buffers &b)
{
	b.ref_bits = std::vector<int>(p.K);
	b.enc_bits = std::vector<int>(p.N);
	b.voltage_level_indexes = std::vector<unsigned>(p.N);
	b.noisy_voltage_levels = std::vector<float>(p.N);
	b.noisy_bits = std::vector<float>(p.N);
	b.dec_bits = std::vector<int>(p.K);
}

void init_utils(const modules &m, utils &u)
{
	// create a sigma noise type
	u.noise = std::unique_ptr<tools::Sigma_asymmetric<>>(new tools::Sigma_asymmetric<>());
	// report the noise values (Es/N0 and Eb/N0)
	u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_noise<>(*u.noise)));
	// report the bit/frame error rates
	u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_BFER<>(*m.monitor)));
	// report the simulation throughputs
	u.reporters.push_back(std::unique_ptr<tools::Reporter>(new tools::Reporter_throughput<>(*m.monitor)));
	// create a terminal that will display the collected data from the reporters
	u.terminal = std::unique_ptr<tools::Terminal_std>(new tools::Terminal_std(u.reporters));
}
