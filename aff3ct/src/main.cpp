#include <iostream>
#include <fstream>
#include <memory>
#include <vector>
#include <string>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct;

struct params
{
	int K = 16384;	// number of information bits
	int N = 18688;			// codeword size

	int fe = 100;			// number of frame errors
	int seed = 0;			// PRNG seed for the AWGN channel
	float ebn0_min = 2.f; // minimum SNR value
	float ebn0_max = 7.9f;	// maximum SNR value
	float ebn0_step = .5f;	// SNR step
	float R;				// code rate (R=K/N)
	float min_sigma = 0.5;	// Set minimum value for individual sigmas (NOT IN USE)
	float skew = 0.45;

	std::vector<float> voltage_levels{-1, 1};
	std::string reader_fpath = "";
	std::string output_fpath = "";
	std::string ldpc_fpath = "";
	int page_type = tools::Flash_reader<float, float>::lower;
	int read_type = tools::Flash_reader<float, float>::soft_single;
	int cell_type = tools::Flash_cell::SLC;

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
	std::unique_ptr<tools::Sigma_asymmetric<>> noise;		 // a sigma noise type
	std::vector<std::unique_ptr<tools::Reporter>> reporters; // list of reporters dispayed in the terminal
	std::unique_ptr<tools::Terminal_std> terminal;			 // manage the output text in the terminal
	std::ofstream output;
};
void init_utils(const params &p, const modules &m, utils &u);
void print_params(std::ostream &stream);

int main(int argc, char **argv)
{
	// get the AFF3CT version
	const std::string v = "v" + std::to_string(tools::version_major()) + "." +
						  std::to_string(tools::version_minor()) + "." +
						  std::to_string(tools::version_release());

	params p;
	init_params(p); // create and initialize the parameters defined by the user
	modules m;
	init_modules(p, m); // create and initialize the modules
	buffers b;
	init_buffers(p, b); // create and initialize the buffers required by the modules
	utils u;
	init_utils(p, m, u); // create and initialize the utils

	print_params(u.output);
	u.terminal->legend(u.output);
	// display the legend in the terminal
	u.terminal->legend();

	// loop over the various SNRs
	std::vector<float> ebn0_range = {7.90643737,  7.69459263,  7.47266019,  7.23966357,  6.99448019, 6.73581115,  6.46214283,  6.17169751,  5.86236912,  5.53163847, 5.17645948,  4.79310398, 4.37694552,  3.9221514 ,  3.42123233, 2.86436372,  2.23832529,  1.52477269,  0.69726958, -0.28414939};
	
	for (auto ebn0_index = 20; ebn0_index >= 0; ebn0_index -= 1)
	{
		// compute the current sigma for the channel noise
		const auto ebn0 = ebn0_range[ebn0_index];
		const auto esn0 = tools::ebn0_to_esn0(ebn0, p.R);
		const auto sigma_ave = tools::esn0_to_sigma(esn0);

		// update the sigma of the modem and the channe
		u.noise->set_sigmas(sigma_ave, ebn0, esn0, p.skew);
		m.channel->set_noise(*u.noise);

		// display the performance (BER and FER) in real time (in a separate thread)
		u.terminal->start_temp_report();

		// run the simulation chain
		while (!m.monitor->fe_limit_achieved() && !u.terminal->is_interrupt())
		{
			m.source->generate(b.ref_bits);
			m.encoder->encode(b.ref_bits, b.enc_bits);
			m.modem->modulate(b.enc_bits, b.voltage_level_indexes);

			m.modem->set_noise(*m.channel);

			m.channel->add_noise(b.voltage_level_indexes.data(), b.noisy_voltage_levels.data());
			m.modem->demodulate(b.noisy_voltage_levels, b.noisy_bits);
			m.decoder->decode_siho(b.noisy_bits, b.dec_bits);
			m.monitor->check_errors(b.dec_bits, b.ref_bits);
		}

		// display the performance (BER and FER) in the terminal
		u.terminal->final_report();
		u.terminal->final_report(u.output);

		// reset the monitor for the next SNR
		m.monitor->reset();
		u.terminal->reset();

		// if user pressed Ctrl+c twice, exit the SNRs loop
		if (u.terminal->is_over())
			break;
	}

	std::cout << "# End of the simulation" << std::endl;
	u.output.close();

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

	const tools::Sparse_matrix H = tools::LDPC_matrix_handler::read(p.ldpc_fpath);
	m.encoder = std::unique_ptr<module::Encoder_LDPC_from_QC<>>(new module::Encoder_LDPC_from_QC<>(p.K, p.N, H));

	tools::Flash_cell cell(p.cell_type);
	tools::Flash_reader<float, float> reader(p.page_type, p.read_type, p.reader_fpath);
	m.modem = std::unique_ptr<module::Modem_flash_page<>>(new module::Modem_flash_page<>(p.N, cell, reader));

	tools::Sigma_asymmetric<float> sigmas;
	sigmas.set_sigmas(1, 2, 2, 0.5);

	m.channel = std::unique_ptr<module::Channel_AWGN_asymmetric<>>(new module::Channel_AWGN_asymmetric<float>(p.N, p.voltage_levels, sigmas));

	m.decoder = std::unique_ptr<module::Decoder_LDPC_BP_flooding_SPA<>>(new module::Decoder_LDPC_BP_flooding_SPA<>(
		p.K,
		p.N,
		50,
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

void init_utils(const params &p, const modules &m, utils &u)
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

	u.output.open(p.output_fpath, std::ios_base::app);
}

void print_params(std::ostream &stream)
{
	// OBS Hardcoded and false parameters
	stream << "[metadata]" << std::endl;
	stream << "ci=on" << std::endl;
	stream << "command=None" << std::endl;
	stream << "title=BCH (7,4,1)" << std::endl;
	stream << "" << std::endl;
	stream << "[trace]" << std::endl;
	stream << "# ----------------------------------------------------" << std::endl;
	stream << "# ---- A FAST FORWARD ERROR CORRECTION TOOLBOX >> ----" << std::endl;
	stream << "# ----------------------------------------------------" << std::endl;
	stream << "# Parameters :" << std::endl;
	stream << "# * Simulation ------------------------------------" << std::endl;
	stream << "#    ** Type                   = BFER" << std::endl;
	stream << "#    ** Type of bits           = int32" << std::endl;
	stream << "#    ** Type of reals          = float32" << std::endl;
	stream << "#    ** Date (UTC)             = 2022-05-10" << std::endl;
	stream << "#    ** Git version            = " << std::endl;
	stream << "#    ** Code type (C)          = QC-LDPC" << std::endl;
	stream << "#    ** Noise range            = 0 -> 10 dB" << std::endl;
	stream << "#    ** Noise type (E)         = EBN0" << std::endl;
	stream << "#    ** Seed                   = 0" << std::endl;
	stream << "#    ** Statistics             = off" << std::endl;
	stream << "#    ** Debug mode             = off" << std::endl;
	stream << "#    ** Multi-threading (t)    = no" << std::endl;
	stream << "#    ** Coset approach (c)     = no" << std::endl;
	stream << "#    ** Coded monitoring       = no" << std::endl;
	stream << "#    ** Bad frames tracking    = off" << std::endl;
	stream << "#    ** Bad frames replay      = off" << std::endl;
	stream << "#    ** Bit rate               = -" << std::endl;
	stream << "#    ** Inter frame level      = -" << std::endl;
	stream << "# * Source ----------------------------------------" << std::endl;
	stream << "#    ** Type                   = RAND" << std::endl;
	stream << "#    ** Implementation         = STD" << std::endl;
	stream << "#    ** Info. bits (K_info)    = 16 384" << std::endl;
	stream << "# * Codec -----------------------------------------" << std::endl;
	stream << "#    ** Type                   = LDPC" << std::endl;
	stream << "#    ** Info. bits (K)         = 16 384" << std::endl;
	stream << "#    ** Codeword size (N_cw)   = 18 688" << std::endl;
	stream << "#    ** Frame size (N)         = 18 688" << std::endl;
	stream << "#    ** Code rate              = 0.877" << std::endl;
	stream << "# * Encoder ---------------------------------------" << std::endl;
	stream << "#    ** Type                   = LDPC" << std::endl;
	stream << "#    ** Systematic             = No" << std::endl;
	stream << "# * Decoder ---------------------------------------" << std::endl;
	stream << "#    ** Type (D)               = BP" << std::endl;
	stream << "#    ** Implementation         = STD" << std::endl;
	stream << "#    ** Systematic             = No" << std::endl;
	stream << "#    ** Galois field order (m) = 2" << std::endl;
	stream << "#    ** Correction power (T)   = 1" << std::endl;
	stream << "# * Modem -----------------------------------------" << std::endl;
	stream << "#    ** Type                   = Flash" << std::endl;
	stream << "#    ** Implementation         = STD" << std::endl;
	stream << "#    ** Bits per symbol        = 1" << std::endl;
	stream << "#    ** Sampling factor        = 1" << std::endl;
	stream << "#    ** Sigma square           = on" << std::endl;
	stream << "# * Channel ---------------------------------------" << std::endl;
	stream << "#    ** Type                   = BiAWGN" << std::endl;
	stream << "#    ** Implementation         = STD" << std::endl;
	stream << "#    ** Complex                = off" << std::endl;
	stream << "#    ** Add users              = off" << std::endl;
	stream << "# * Monitor ---------------------------------------" << std::endl;
	stream << "#    ** Frame error count (e)  = 100" << std::endl;
	stream << "#    ** Compute mutual info    = no" << std::endl;
	stream << "# * Terminal --------------------------------------" << std::endl;
	stream << "#    ** Enabled                = yes" << std::endl;
	stream << "#    ** Frequency (ms)         = 500" << std::endl;
	stream << "#" << std::endl;
	stream << "# The simulation is running..." << std::endl;
}