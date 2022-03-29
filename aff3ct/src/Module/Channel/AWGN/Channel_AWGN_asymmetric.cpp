#include <algorithm>
#include <string>
#include <vector>
#include <map>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct;
using namespace aff3ct::module;
using std::vector;

template <typename R>
Channel_AWGN_asymmetric<R>::Channel_AWGN_asymmetric(const int N, std::unique_ptr<tools::Gaussian_gen<R>> &&_ng, const bool add_users,
							  const tools::Sigma<R> &total_noise, const int n_frames)
	: Channel<R>(N, total_noise, n_frames),
	  add_users(add_users),
	  noise_generator(std::move(_ng))
{
	const std::string name = "Channel_AWGN_asymmetric";
	this->set_name(name);

	if (this->noise_generator == nullptr)
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, "'noise_generator' can't be NULL.");
	
	this->noise_map = new tools::Sigma_asymmetric<R>(total_noise);
}

template <typename R>
Channel_AWGN_asymmetric<R>::Channel_AWGN_asymmetric(const int N, const int seed, const bool add_users, const tools::Sigma<R> &total_noise, const int n_frames)
	: Channel_AWGN_asymmetric<R>(N, std::unique_ptr<tools::Gaussian_noise_generator_std<R>>(new tools::Gaussian_noise_generator_std<R>(seed)),
					  add_users, total_noise, n_frames)
{
}

template <typename R>
void Channel_AWGN_asymmetric<R>::add_noise(const R *X_N, R *Y_N, const int frame_id)
{
	this->check_noise();

	if (add_users) //&& this->n_frames > 1)
	{
		// Nothing happens when add_user is true
		if (frame_id != -1)
		{
			std::stringstream message;
			message << "'frame_id' has to be equal to -1 ('frame_id' = " << frame_id << ").";
			throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
		}
	}
	else
	{
		const auto f_start = (frame_id < 0) ? 0 : frame_id % this->n_frames;
		const auto f_stop = (frame_id < 0) ? this->n_frames : f_start + 1;

		if (frame_id < 0)
			noise_generator->generate(this->noise, this->n->get_noise(level));
		else
			// Not correctly implemented
			noise_generator->generate(this->noise.data() + f_start * this->N, this->N, this->n->get_noise());

		for (auto f = f_start; f < f_stop; f++) {
			for (auto n = 0; n < this->N; n++) {
				noise_generator->generate(&this->noise[f * this->N + n],1,sigma[Y_N], this->n->get_noise());
				Y_N[f * this->N + n] = X_N[f * this->N + n] + this->noise[f * this->N + n];
			}

		}
			
				

		/*
		for (auto f = 0; f < 1; f++) {
		for (size_t n = 0; n < noise.size(); n++) {
		    voltage_adress_map[c[n]].push_back(&noise[n]);
		}
		// For each voltage level, apply a distribution
		auto map_iterator = voltage_adress_map.begin();
		for (size_t i = 0; i < voltage_adress_map.size(); i++) {
		    generate(map_iterator->second[0], map_iterator->second.size(), sigmas[i], 0.0f);
		    map_iterator++;
		}
		
		// Add the noise to the modulated data
		for (size_t n = 0; n < noise.size(); n++) {
		    c[n] += noise[n];
		    cout << c[n] << endl;
		}
		}
		*/
	}
}

template <typename R>
void Channel_AWGN_asymmetric<R>::check_noise()
{
	Channel<R>::check_noise();

	this->n->is_of_type_throw(tools::Noise_type::SIGMA);
}
// ==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Channel_AWGN_asymmetric<R_32>;
template class aff3ct::module::Channel_AWGN_asymmetric<R_64>;
#else
template class aff3ct::module::Channel_AWGN_asymmetric<R>;
#endif
// ==================================================================================== explicit template instantiation
