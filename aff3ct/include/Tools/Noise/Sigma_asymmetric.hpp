/*!
 * \file
 * \brief Class tools::Sigma_asymmetric.
 */
#ifndef SIGMA_ASYMMETRIC_HPP_
#define SIGMA_ASYMMETRIC_HPP_

#include <utility>
#include <map>

#include "Tools/Noise/Noise.hpp"

namespace aff3ct
{
	namespace tools
	{

		template <typename R = float>
		class Sigma_asymmetric : public Sigma<R>
		{
		public:
			Sigma_asymmetric();
			explicit Sigma_asymmetric(R tot_noise, R min_sigma);
			Sigma_asymmetric(R tot_noise, R min_sigma, R tot_ebn0, R tot_esn0);
			template <typename T>
			explicit Sigma_asymmetric(const Sigma_asymmetric<T> &other);
			virtual ~Sigma_asymmetric() = default;

			R get_sigma(R voltage_level) const;
			std::map<R, R> get_sigmas() const;
			void generate_sigmas();

			virtual Noise_type get_type() const;

			bool has_sigmas() const noexcept;

			virtual void copy(const Sigma_asymmetric &other); // set this noise as the 'other' one

			virtual Sigma_asymmetric<R> *clone() const;

		protected:
			std::map<R, R> sigma_map; // Maps a sigma to a given voltage level
			bool sigmas_exist = false;
			R min_sigma;
			tools::Random_sigma_generator<R> sigma_generator;

			virtual void check();
		};

	}
}

#endif // SIGMA_ASYMMETRIC_HPP_