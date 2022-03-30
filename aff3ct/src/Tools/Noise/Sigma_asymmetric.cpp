#include <sstream>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct;
using namespace aff3ct::tools;

template<typename R>
Sigma_asymmetric<R>::Sigma_asymmetric()
: Sigma<R>()
{

}

template<typename R>
Sigma_asymmetric<R>::Sigma_asymmetric(R tot_noise)
: Sigma<R>(tot_noise)
{
	sigma_map = std::map<R,R>();
	this -> generate_sigmas();
}

template<typename R>
Sigma_asymmetric<R>::Sigma_asymmetric(R tot_noise, R tot_ebn0, R tot_esn0)
: Sigma<R>(tot_noise, tot_ebn0, tot_esn0)
{
	sigma_map = std::map<R,R>();
	this -> generate_sigmas();
	sigmas_exist = true;
}

template <typename R>
template <typename T>
Sigma_asymmetric<R>::
Sigma_asymmetric(const Sigma_asymmetric<T>& other)
: Sigma<R>(other)
{
	if (other.has_sigmas())
	{
		sigma_map = std::map<R,R>();
		auto other_sigma_map = other.get_sigmas();
		for (auto pair: other_sigma_map) {
			R voltage_level = (R) pair.first;
			R sigma = (R) pair.second;
			sigma_map[voltage_level] = sigma; 
		}
	}
	else
	{
		sigma_map = std::map<R,R>();
		sigmas_exist = false;
	}
}

template <typename R>
void Sigma_asymmetric<R>::
copy(const Sigma_asymmetric<R>& other)
{
	sigma_map = other.sigma_map;
	sigmas_exist = other.sigma_exist;
	Sigma<R>::copy(other);
}


template<typename R>
void Sigma_asymmetric<R>::generate_sigmas()
{
	return _ebn0.second;
}

template<typename R>
R Sigma_asymmetric<R>::get_sigma(R voltage_level) const {
	return sigma_map[voltage_level];
}

template<typename R>
std::map<R,R> Sigma_asymmetric<R>::get_sigmas() const {
	return sigma_map;
}


template<typename R>
bool Sigma_asymmetric<R>::has_sigmas() const noexcept
{
	return sigmas_exist;
}

template <typename R>
void Sigma_asymmetric<R>::
check()
{
	auto n = this->get_noise();
	if (n <= (R)0)
	{
		std::stringstream message;
		message << "The SIGMA_ASYMMETRIC total noise '_n' has to be greater than 0 ('_n' = " << n << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}
}

template<typename R>
Noise_type Sigma_asymmetric<R>::get_type() const
{
	// Add SIGMA_ASYMMETRIC to noise types
	return Noise_type::SIGMA;
}

template<typename R>
Sigma_asymmetric<R>* Sigma_asymmetric<R>::clone() const
{
	return new Sigma_asymmetric<R>(*this);
}



// ==================================================================================== explicit template instantiation
template class aff3ct::tools::Sigma_asymmetric<float>;
template class aff3ct::tools::Sigma_asymmetric<double>;

template aff3ct::tools::Sigma_asymmetric<double>::Sigma_asymmetric(const Sigma_asymmetric<float >&);
template aff3ct::tools::Sigma_asymmetric<float >::Sigma_asymmetric(const Sigma_asymmetric<double>&);
// ==================================================================================== explicit template instantiation
