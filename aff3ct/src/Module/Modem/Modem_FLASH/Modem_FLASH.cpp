#include <string>
#include <algorithm>
#include <type_traits>

#include <aff3ct.hpp> 
#include <aff3ct_extension.hpp>

using namespace aff3ct;
using namespace aff3ct::module;

template <typename B, typename R, typename Q>
Modem_FLASH<B,R,Q>
::Modem_FLASH(const int N, const tools::Noise<R>& noise, const int n_frames)
:Modem<B,R,Q>(N, noise, n_frames) 
{

}

template <typename B, typename R, typename Q>
void Modem_FLASH<B,R,Q>
::set_noise(const tools::Noise<R>& noise)
{
	Modem<B,R,Q>::set_noise(noise);

	this->n->is_of_type_throw(tools::Noise_type::SIGMA);
}


template <typename B, typename R, typename Q>
void Modem_FLASH<B,R,Q>
::_modulate(const B *X_N1, R *X_N2, const int frame_id)
{
    //TODO: Replace dummy function
	std::copy(X_N1, X_N1, X_N2);
}

template <typename B,typename R, typename Q>
void Modem_FLASH<B,R,Q>
::_filter(const R *Y_N1, R *Y_N2, const int frame_id)
{
	std::copy(Y_N1, Y_N1 + this->N_fil, Y_N2);
}

template <typename B, typename R, typename Q>
void Modem_FLASH<B,R,Q>
::_demodulate(const Q *Y_N1, Q *Y_N2, const int frame_id)
{
    //TODO: Replace dummy function
	std::copy(Y_N1,Y_N1, Y_N2);
}



// ==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::module::Modem_FLASH<B_8,R_8,R_8>;
template class aff3ct::module::Modem_FLASH<B_8,R_8,Q_8>;
template class aff3ct::module::Modem_FLASH<B_16,R_16,R_16>;
template class aff3ct::module::Modem_FLASH<B_16,R_16,Q_16>;
template class aff3ct::module::Modem_FLASH<B_32,R_32,R_32>;
template class aff3ct::module::Modem_FLASH<B_64,R_64,R_64>;
#else
template class aff3ct::module::Modem_FLASH<B,R,Q>;
#if !defined(AFF3CT_32BIT_PREC) && !defined(AFF3CT_64BIT_PREC)
template class aff3ct::module::Modem_FLASH<B,R,R>;
#endif
#endif
// ==================================================================================== explicit template instantiation
