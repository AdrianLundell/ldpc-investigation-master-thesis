#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>
#include <math.h>

using namespace aff3ct;
using namespace aff3ct::tools;

template <typename R>
Sigma_asymmetric<R>::Sigma_asymmetric()
	: Sigma<R>()
{
	this->sigmas = std::vector<R>(2, (R)0.0);
}

template <typename R>
Sigma_asymmetric<R>::Sigma_asymmetric(const Sigma_asymmetric<R> &other)
	: Sigma_asymmetric<R>()
{
	this->sigmas = other.sigmas;
	this->skew = other.skew;
	this->sigmas_exist = other.sigmas_exist;
	this->set_noise(other.get_noise(), other.get_ebn0(), other.get_esn0());
}

template <typename R>
void Sigma_asymmetric<R>::set_sigmas(R sigma_ave, R ebn0, R esn0, R skew)
{
	this->set_noise(sigma_ave, ebn0, esn0);
	this->skew = skew;

	this->sigmas[0] = (2 * sigma_ave) * (1 - skew);
	this->sigmas[1] = (2 * sigma_ave) * skew;

	this->sigmas_exist = true;
	this->check();
}

template <typename R>
R Sigma_asymmetric<R>::get_sigma(unsigned voltage_index) const
{
	return sigmas[voltage_index];
}

template <typename R>
R Sigma_asymmetric<R>::get_sigma_ave() const
{
	return this->get_noise();
}

template <typename R>
std::vector<R> Sigma_asymmetric<R>::get_sigmas() const
{
	return sigmas;
}

template <typename R>
bool Sigma_asymmetric<R>::has_sigmas() const noexcept
{
	return sigmas_exist;
}

template <typename R>
R Sigma_asymmetric<R>::get_threshold_noise(unsigned threshold_index) const
{
	return (R)-10 * log10(pow(sigmas[threshold_index], (R)2.0) + pow(sigmas[threshold_index + 1], (R)2.0));
}

template <typename R>
void Sigma_asymmetric<R>::check() const
{
	R n = this->get_noise();

	if (skew <= (R)0 || skew >= (R)1)
	{
		std::stringstream message;
		message << "The SIGMA_ASYMMETRIC skewness 'skew' has to be greater than 0 and smaller than 1" << skew << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}
}

template <typename R>
Noise_type Sigma_asymmetric<R>::get_type() const
{
	// Add SIGMA_ASYMMETRIC to noise types
	return Noise_type::SIGMA;
}

/*
template <typename R>
Sigma_asymmetric<R> *Sigma_asymmetric<R>::clone() const
{
	return new Sigma_asymmetric<R>(*this);
}
*/

// ==================================================================================== explicit template instantiation
template class aff3ct::tools::Sigma_asymmetric<float>;
template class aff3ct::tools::Sigma_asymmetric<double>;

// template aff3ct::tools::Sigma_asymmetric<double>::Sigma_asymmetric(const Sigma_asymmetric<float> &);
// template aff3ct::tools::Sigma_asymmetric<float>::Sigma_asymmetric(const Sigma_asymmetric<double> &);
//  ==================================================================================== explicit template instantiation
