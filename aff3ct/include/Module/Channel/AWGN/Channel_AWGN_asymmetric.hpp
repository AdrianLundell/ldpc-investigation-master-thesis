/*!
 * \file
 * \brief Class module::Channel_AWGN_asymmetric.
 */
#ifndef CHANNEL_AWGN_ASYMMETRIC_HPP_
#define CHANNEL_AWGN_ASYMMETRIC_HPP_

#include <memory>
#include <vector>
#include <aff3ct.hpp>

namespace aff3ct
{
	namespace module
	{
		template <typename R = float>
		class Channel_AWGN_asymmetric : public Channel<R>
		{
		private:
			const bool add_users;
			std::unique_ptr<tools::Gaussian_gen<R>> noise_generator;
			std::unique_ptr<tools::Sigma_asymmetric<R>> noise_map;

		public:
			Channel_AWGN_asymmetric(const int N, std::unique_ptr<tools::Gaussian_gen<R>> &&noise_generator,
						 const bool add_users = false,
						 const tools::Sigma<R> &total_noise = tools::Sigma<R>((R)1),
						 const int n_frames = 1);

			explicit Channel_AWGN_asymmetric(const int N, const int seed = 0, const bool add_users = false,
								  const tools::Sigma<R> &total_noise = tools::Sigma<R>((R)1),
								  const int n_frames = 1);

			virtual ~Channel_AWGN_asymmetric() = default;

			void add_noise(const R *X_N, R *Y_N, const int frame_id = -1);
			using Channel<R>::add_noise;

		protected:
			virtual void check_noise();
		};
	}
}

#endif /* CHANNEL_AWGN_asymmetric_HPP_ */
