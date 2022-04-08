/*!
 * \file
 * \brief Class tools::Random_sigma_generator.
 */
#ifndef RANDOM_SIGMA_GENERATOR_HPP_
#define RANDOM_SIGMA_GENERATOR_HPP_

#include <vector>
#include <memory>
#include <random>

#include "Tools/Algo/Draw_generator/Draw_generator.hpp"

namespace aff3ct
{
	namespace tools
	{
		template <typename R = float>
		class Random_sigma_generator : public Draw_generator<R>
		{
		private:
			//	std::minstd_rand            rd_engine; // LCG
			std::mt19937 rd_engine; // Mersenne Twister 19937
			std::uniform_real_distribution<R> uniform_real_dist;

		public:
			Random_sigma_generator() = default;
			explicit Random_sigma_generator(const int seed = 0);

			virtual ~Random_sigma_generator() = default;

			virtual void generate(std::vector<R> &sigmas, const R sigma_tot, const R min_sigma);

			virtual void set_seed(const int seed);
		};

	}
}

#endif /* RANDOM_SIGMA_GENERATOR_HPP_ */
