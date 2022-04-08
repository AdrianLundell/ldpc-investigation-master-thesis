/*!
 * \file
 * \brief Class module::Channel_AWGN_asymmetric.
 */
#ifndef CHANNEL_AWGN_ASYMMETRIC_HPP_
#define CHANNEL_AWGN_ASYMMETRIC_HPP_

#include <memory>
#include <vector>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

namespace aff3ct
{
	namespace module
	{
		template <typename R = float>
		class Channel_AWGN_asymmetric : public Channel<R>
		{
		private:
			//std::unique_ptr<tools::Gaussian_noise_generator_std<R>> noise_generator;
			tools::Gaussian_noise_generator_std<R> noise_generator;
			tools::Sigma_asymmetric<R> sigmas;
			std::vector<R> voltage_levels;

		public:
			Channel_AWGN_asymmetric(const int N,
						std::vector<R> &voltage_levels,
						 const int seed = 0,
						 const int n_frames = 1);

			virtual ~Channel_AWGN_asymmetric() = default;

			virtual void set_noise(tools::Sigma_asymmetric<R>& noise);

			virtual void add_noise(const unsigned *voltage_level_indexes, R *noisy_voltage_levels, const int frame_id = -1);
			//using Channel<R>::add_noise;
			R get_snr(const unsigned threshold_index) const;
			R get_sigma_ratio(const unsigned threshold_index) const;

		protected:
			virtual void check_noise();
		};
	}
}

#endif /* CHANNEL_AWGN_asymmetric_HPP_ */
