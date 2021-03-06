#include <algorithm>
#include <string>
#include <vector>
#include <memory>
#include <math.h>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct;
using namespace aff3ct::module;
using std::vector;

template <typename R>
Channel_AWGN_asymmetric<R>::Channel_AWGN_asymmetric(const int N,
													std::vector<R> v_l,
													const tools::Sigma_asymmetric<R> &_s,
													const int seed,
													const int n_frames)
	: Channel<R>(N, n_frames),

	  voltage_levels(v_l),
	  sigmas(_s)
{
	const std::string name = "Channel_AWGN_asymmetric";
	this->set_name(name);
	// noise_generator = std::unique_ptr<tools::Gaussian_noise_generator_std<R>>(new tools::Gaussian_noise_generator_std<R>(seed));
	noise_generator = tools::Gaussian_noise_generator_std<R>(seed);
}

template <typename R>
Channel_AWGN_asymmetric<R>::Channel_AWGN_asymmetric(const int N,
													std::vector<R> v_l,
													const int seed,
													const int n_frames)
	: Channel_AWGN_asymmetric(N, v_l, tools::Sigma_asymmetric<R>(), seed, n_frames)
{
}

template <typename R>
void Channel_AWGN_asymmetric<R>::set_noise(tools::Sigma_asymmetric<R> &_s)
{
	this->sigmas = _s;
}

template <typename R>
R Channel_AWGN_asymmetric<R>::get_sigma(unsigned voltage_level_index) const
{
	return sigmas.get_sigma(voltage_level_index);
}

template <typename R>
R Channel_AWGN_asymmetric<R>::get_sigma_ave() const
{
	return sigmas.get_sigma_ave();
}

template <typename R>
void Channel_AWGN_asymmetric<R>::add_noise(const unsigned *voltage_level_indexes, R *noisy_voltage_levels, const int frame_id)
{
	const auto f_start = (frame_id < 0) ? 0 : frame_id % this->n_frames;
	const auto f_stop = (frame_id < 0) ? this->n_frames : f_start + 1;

	for (auto f = f_start; f < f_stop; f++)
	{
		for (auto n = 0; n < this->N; n++)
		{
			unsigned current_idx = f * this->N + n;
			unsigned voltage_level_index = voltage_level_indexes[current_idx];
			// Generate noise distributed around the current voltage level
			noise_generator.generate(noisy_voltage_levels, (unsigned)1,
									 this->sigmas.get_sigma(voltage_level_index),
									 this->voltage_levels[voltage_level_index]);

			noisy_voltage_levels++;
		}
	}
}

template <typename R>
R Channel_AWGN_asymmetric<R>::get_snr(const unsigned threshold_index) const
{
	// R signal_term = 10 * log10(pow(voltage_levels[threshold_index] - voltage_levels[threshold_index + 1], (R)2.0));
	// R noise_term = sigmas.get_threshold_noise(threshold_index);
	// R result = signal_term + noise_term;

	return sigmas.get_ebn0();
}

template <typename R>
R Channel_AWGN_asymmetric<R>::get_skew(const unsigned threshold_index) const
{
	if (threshold_index > 0)
	{
		std::stringstream message;
		message << "The SIGMA_ASYMMETRIC only contains two adjacent distributions and thus, the threshold index cannot be above 0. It is: " << threshold_index << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	R sigma1 = sigmas.get_sigma(0);
	R sigma2 = sigmas.get_sigma(1);
	return (sigma2) / (sigma1 + sigma2);
}

template <typename R>
void Channel_AWGN_asymmetric<R>::check_noise()
{
	Channel<R>::check_noise();
	this->n->is_of_type_throw(tools::Noise_type::SIGMA);
}

//==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Channel_AWGN_asymmetric<R_32>;
template class aff3ct::module::Channel_AWGN_asymmetric<R_64>;
#else
template class aff3ct::module::Channel_AWGN_asymmetric<R>;
#endif
//==================================================================================== explicit template instantiation
