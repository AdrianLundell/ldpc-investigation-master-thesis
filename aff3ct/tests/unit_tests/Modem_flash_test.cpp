// #include <gtest/gtest.h>

// #include <memory>
// #include <vector>
// #include <string>

// #include <aff3ct.hpp>
// #include <aff3ct_extension.hpp>

// using namespace aff3ct;

// struct params
// {
// 	int K = 8;				 // number of information bits
// 	int N = 8;			     // codeword size
// 	int fe = 100;			 // number of frame errors
// 	int seed = 0;			 // PRNG seed for the AWGN channel
// 	float ebn0_min = 0.00f;	 // minimum SNR value
// 	float ebn0_max = 0.01f;  // maximum SNR value
// 	float ebn0_step = 1.00f; // SNR step
// 	float R;				 // code rate (R=K/N)
// };
// void init_params(params &p);

// struct modules
// {
// 	std::unique_ptr<module::Source_random<>> source;
// 	std::unique_ptr<module::Encoder_repetition_sys<>> encoder;
// 	std::unique_ptr<module::Modem_flash<>> modem;
// 	std::unique_ptr<module::Channel_Test<>> channel;
// 	std::unique_ptr<module::Decoder_repetition_std<>> decoder;
// 	std::unique_ptr<module::Monitor_BFER<>> monitor;
// };
// void init_modules(const params &p, modules &m);

// struct buffers
// {
// 	std::vector<int> ref_bits;
// 	std::vector<int> enc_bits;
// 	std::vector<float> symbols;
// 	std::vector<float> noisy_symbols;
// 	std::vector<float> LLRs;
// 	std::vector<int> dec_bits;
// };
// void init_buffers(const params &p, buffers &b);

// struct utils
// {
// 	std::unique_ptr<tools::Sigma<>> noise;
// };
// void init_utils(const modules &m, utils &u);

// void init_params(params &p)
// {
// 	p.R = (float)p.K / (float)p.N;
// }

// void init_modules(const params &p, modules &m)
// {
// 	m.source = std::unique_ptr<module::Source_random<>>(new module::Source_random<>(p.K));
// 	m.encoder = std::unique_ptr<module::Encoder_repetition_sys<>>(new module::Encoder_repetition_sys<>(p.K, p.N));
// 	m.modem = std::unique_ptr<module::Modem_flash<>>(new module::Modem_flash<>(p.N, 
// 			  std::unique_ptr<tools::Constellation<float>>(new tools::Constellation_flash<float>("SLC_bpsk_voltage_levels.txt")),
// 			  std::unique_ptr<tools::Thresholder<float>>(new tools::Thresholder_soft<float>("test_thresholds.txt"))));
// 	m.channel = std::unique_ptr<module::Channel_Test<>>(new module::Channel_Test<>(p.N, p.seed));
// 	m.decoder = std::unique_ptr<module::Decoder_repetition_std<>>(new module::Decoder_repetition_std<>(p.K, p.N));
// 	m.monitor = std::unique_ptr<module::Monitor_BFER<>>(new module::Monitor_BFER<>(p.K, p.fe));
// };

// void init_buffers(const params &p, buffers &b)
// {
// 	b.ref_bits = std::vector<int>(p.K);
// 	b.enc_bits = std::vector<int>(p.N);
// 	b.symbols = std::vector<float>(p.N);
// 	b.noisy_symbols = std::vector<float>(p.N);
// 	b.LLRs = std::vector<float>(p.N);
// 	b.dec_bits = std::vector<int>(p.K);
// }

// void init_utils(const modules &m, utils &u)
// {
// 	// create a sigma noise type
// 	u.noise = std::unique_ptr<tools::Sigma<>>(new tools::Sigma<>());
// }

// class ModemFlashTest : public ::testing::Test {
// 	protected:
// 		void SetUp() override {
// 			params p;
// 			init_params(p); // create and initialize the parameters defined by the user
// 			modules m;
// 			init_modules(p, m); // create and initialize the modules
// 			buffers b;
// 			init_buffers(p, b); // create and initialize the buffers required by the modules
// 			utils u;
// 			init_utils(m, u); // create and initialize the utils

// 			m.source->generate(b.ref_bits);
// 			m.encoder->encode(b.ref_bits, b.enc_bits);
// 		}
// }

// TEST_F(ModemFlashTest, ModulationTests) {
// 	m.modem->modulate(b.enc_bits, b.symbols);
// 		m.channel->add_noise(b.symbols, b.noisy_symbols);
// 		m.modem->demodulate(b.noisy_symbols, b.LLRs);
// 		ASSERT()
// }