#include <sstream>
#include <math.h>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct;
using namespace aff3ct::tools;

template <typename R>
Sigma_asymmetric<R>::Sigma_asymmetric()
	: Sigma<R>()
{
	this->sigma_generator = tools::Random_sigma_generator<R>();
}

template <typename R>
Sigma_asymmetric<R>::Sigma_asymmetric(const Sigma_asymmetric<R> &other)
	: Sigma_asymmetric<R>()
{
	this->min_sigma = other.min_sigma;
	this->sigmas = other.sigmas;
	this->sigma_ratios = other.sigma_ratios;
	this->set_noise(other.get_noise(), other.get_ebn0(), other.get_esn0());
}

template <typename R>
void Sigma_asymmetric<R>::set_sigmas(R sigma_tot, unsigned n_sigmas, R m_s, R ebn0, R esn0)
{
	this->min_sigma = m_s;
	this->sigmas = std::vector<R>(n_sigmas, (R)0.0);
	this->sigma_ratios = std::vector<R>(n_sigmas - 1, (R)0.0);
	this->set_noise(sigma_tot, ebn0, esn0);
	this->sigmas_exist = true;
	this->check();
}

template <typename R>
void Sigma_asymmetric<R>::set_sigmas(std::vector<R>& sigmas)
{
	R min = (R) sigmas[0];
	R tot = (R) 0;
	R ebn0, esn0;

	for (auto sigma : sigmas)
	{
		tot += sigma;
		if (sigma < min) min = sigma;
	}

	this->min_sigma = min;
	this->sigmas = sigmas;
	this->sigma_ratios = std::vector<R>(sigmas.size() - 1, (R)0.0);
	this->set_noise(tot, ebn0, esn0);
	this->sigmas_exist = true;
	this->check();

	this->compute_sigma_ratios();
}

template <typename R>
R Sigma_asymmetric<R>::get_sigma(unsigned voltage_index) const
{
	return sigmas[voltage_index];
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
void Sigma_asymmetric<R>::generate_sigmas()
{
	if (this->sigmas_exist)
	{
		this->sigma_generator.generate(sigmas, this->get_noise(), min_sigma);
		this->compute_sigma_ratios();
	}
	else
	{
		std::stringstream message;
		message << "The SIGMA_ASYMMETRIC parameters have not been set. Call set_sigmas() before calling this method.";
		throw tools::logic_error(__FILE__, __LINE__, __func__, message.str());
	}
}

template <typename R>
void Sigma_asymmetric<R>::compute_sigma_ratios()
{
	for (unsigned i = 0; i < sigma_ratios.size(); i++)
	{
		this->sigma_ratios[i] = pow(sigmas[i], (R)2.0) / pow(sigmas[i + 1], (R)2.0);
	}
}

template <typename R>
void Sigma_asymmetric<R>::set_seed(const int seed)
{
	this->sigma_generator.set_seed(seed);
}

template <typename R>
R Sigma_asymmetric<R>::get_threshold_noise(unsigned threshold_index) const
{
	return (R)-10 * log10(pow(sigmas[threshold_index], (R)2.0) + pow(sigmas[threshold_index + 1], (R)2.0));
}

template <typename R>
R Sigma_asymmetric<R>::get_ratio(unsigned threshold_index) const
{
	R sigma_i_sqrd = pow(sigmas[threshold_index], (R)2.0);
	R sigma_ii_sqrd = pow(sigmas[threshold_index + 1], (R)2.0);
	return sigma_i_sqrd / sigma_ii_sqrd + sigma_ii_sqrd;
}

template <typename R>
void Sigma_asymmetric<R>::check() const
{
	R n = this->get_noise();
	if (n <= (R)0)
	{
		std::stringstream message;
		message << "The SIGMA_ASYMMETRIC total sigma '_n' has to be greater than 0 ('_n' = " << n << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	if (min_sigma <= (R)0)
	{
		std::stringstream message;
		message << "The SIGMA_ASYMMETRIC minimum sigma 'min_sigma' has to be greater than 0 ('min_sigma' = " << min_sigma << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}

	unsigned n_sigmas = sigmas.size();
	R min_sigma_tot = min_sigma * n_sigmas;
	if (n < min_sigma_tot)
	{
		std::stringstream message;
		message << "The SIGMA_ASYMMETRIC total sigma '_n' has to be greater than " << n_sigmas << "*min_sigma  ('n' = " << n << ", 'min_sigma' =" << min_sigma << ").";
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
