/*!
 * \file
 * \brief Class module::Channel_Test.
 */
#ifndef CHANNEL_TEST_HPP_
#define CHANNEL_TEST_HPP_

#include <memory>
#include <aff3ct.hpp>

namespace aff3ct
{
	namespace module
	{
		template <typename R = float>
		class Channel_Test : public Channel<R>
		{
		private:
			const bool add_users;
			std::unique_ptr<tools::Gaussian_gen<R>> noise_generator;

		public:
			Channel_Test(const int N, std::unique_ptr<tools::Gaussian_gen<R>> &&noise_generator,
						 const bool add_users = false,
						 const tools::Sigma<R> &noise = tools::Sigma<R>(),
						 const int n_frames = 1);

			explicit Channel_Test(const int N, const int seed = 0, const bool add_users = false,
								  const tools::Sigma<R> &noise = tools::Sigma<R>(),
								  const int n_frames = 1);

			virtual ~Channel_Test() = default;

			void add_noise(const R *X_N, R *Y_N, const int frame_id = -1);
			using Channel<R>::add_noise;

		protected:
			virtual void check_noise();
		};
	}
}

#endif /* CHANNEL_Test_HPP_ */
