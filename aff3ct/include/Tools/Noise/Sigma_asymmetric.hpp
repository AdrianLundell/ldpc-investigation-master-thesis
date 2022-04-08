/*!
 * \file
 * \brief Class tools::Sigma_asymmetric.
 */
#ifndef SIGMA_ASYMMETRIC_HPP_
#define SIGMA_ASYMMETRIC_HPP_

#include <utility>
#include <vector>

#include "Tools/Noise/Noise.hpp"
#include <aff3ct_extension.hpp>

namespace aff3ct
{
	namespace tools
	{

		template <typename R = float>
		class Sigma_asymmetric : public Sigma<R>
		{
		public:
			Sigma_asymmetric();
			Sigma_asymmetric(Sigma_asymmetric<R> &other);
			virtual ~Sigma_asymmetric() = default;
			void set_sigmas(R sigma_tot, unsigned n_sigmas, R sigma_min, R ebn0, R esn0);
			R get_sigma(unsigned voltage_level_index) const;
			std::vector<R> get_sigmas() const;
			bool has_sigmas() const noexcept;
			void generate_sigmas();
			R get_threshold_noise(unsigned threshold_index) const;
			R get_ratio(unsigned threshold_index) const;
			
			virtual Noise_type get_type() const;
			
			void check() const;
			virtual Sigma_asymmetric<R> *clone() const;

		protected:
			void compute_sigma_ratios();
			std::vector<R> sigmas; // Maps a sigma to a given voltage level
			std::vector<R> sigma_ratios;
			bool sigmas_exist = false;
			R min_sigma;
			tools::Random_sigma_generator<R> sigma_generator;

			virtual void check();
		};

	}
}

#endif // SIGMA_ASYMMETRIC_HPP_