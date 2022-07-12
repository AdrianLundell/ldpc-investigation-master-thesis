/*!
 * \file
 * \brief Class tools::Sigma_asymmetric.
 *
 * \section Definition:
 * Sigma_asymmetric represents two adjacent normal distributions that are asymmetric. The class inherits from Sigma.
 * sigma_ave, ebn0, esn0 represent the channel parameters if the distributions were symmetric.
 * sigma for each distribution is computed as (2*sigma_ave)*skew for the leftmost distribution and (2*sigma_ave)*(1-skew) for the rightmost distribution.
 */
#ifndef SIGMA_ASYMMETRIC_HPP_
#define SIGMA_ASYMMETRIC_HPP_

#include <utility>
#include <vector>

#include <aff3ct.hpp>
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
			explicit Sigma_asymmetric(const Sigma_asymmetric<R> &other);
			virtual ~Sigma_asymmetric() = default;
			void set_sigmas(R sigma_ave, R ebn0, R esn0, R skew);
			void set_skewness(std::vector<R> &sigmas);
			R get_sigma(unsigned voltage_level_index) const;
			R get_sigma_ave() const;
			std::vector<R> get_sigmas() const;
			bool has_sigmas() const noexcept;
			R get_threshold_noise(unsigned threshold_index) const;
			R get_ratio(unsigned threshold_index) const;

			virtual Noise_type get_type() const;

			virtual void check() const;
			// virtual Sigma_asymmetric<R> *clone() const;

		protected:
			std::vector<R> sigmas; // Maps a sigma to a given voltage level
			R skew;
			bool sigmas_exist = false;
		};

	}
}

#endif // SIGMA_ASYMMETRIC_HPP_