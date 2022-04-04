#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <algorithm>
#include <limits>
#include <string>
#include <sstream>
#include <type_traits>
#include <complex>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

namespace aff3ct
{
namespace module
{
template <typename B, typename R, typename Q, tools::proto_max<Q> MAX>
Modem_flash<B,R,Q,MAX>
::Modem_flash(const int N, std::unique_ptr<const tools::Constellation<R>>&& _cstl, 
				std::unique_ptr<const tools::Thresholder<R>>&& _thresholder,
				const tools::Noise<R>& noise,
				const int n_frames)
: Modem<B,R,Q>(N,
              (int)(std::ceil((float)N / (float)_cstl->get_n_bits_per_symbol())), // N_mod
              noise,
              n_frames),
  cstl           (std::move(_cstl)),
  thresholder    (std::move(_thresholder)),
  bits_per_symbol(cstl->get_n_bits_per_symbol()),
  nbr_symbols    (cstl->get_n_symbols())
{
	const std::string name = "Modem_flash<" + cstl->get_name() + ">";
	this->set_name(name);

	if (cstl == nullptr)
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, "No constellation given ('cstl' = nullptr).");


	if (thresholder == nullptr)
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, "No thresholder given ('thresholder' = nullptr).");
}

template <typename B, typename R, typename Q, tools::proto_max<Q> MAX>
void Modem_flash<B,R,Q,MAX>
::set_noise(const tools::Noise<R>& noise)
{
	Modem<B,R,Q>::set_noise(noise);

	this->n->is_of_type_throw(tools::Noise_type::SIGMA);
}

template <typename B, typename R, typename Q, tools::proto_max<Q> MAX>
int Modem_flash<B,R,Q,MAX>
::size_mod(const int N, const tools::Constellation<R>& c)
{
	return Modem<B,R,Q>::get_buffer_size_after_modulation(N, c.get_n_bits_per_symbol(), 0, 1, false);
}

template <typename B, typename R, typename Q, tools::proto_max<Q> MAX>
int Modem_flash<B,R,Q,MAX>
::size_fil(const int N, const tools::Constellation<R>& c)
{
	return Modem<B,R,Q>::get_buffer_size_after_filtering(N, c.get_n_bits_per_symbol(), 0, 1, false);
}

template <typename B,typename R, typename Q, tools::proto_max<Q> MAX>
void Modem_flash<B,R,Q,MAX>
::_filter(const R *Y_N1, R *Y_N2, const int frame_id)
{
	std::copy(Y_N1, Y_N1 + this->N_fil, Y_N2);
}

template <typename B,typename R, typename Q, tools::proto_max<Q> MAX>
void Modem_flash<B,R,Q,MAX>
::_modulate(const B *X_N1, R *X_N2, const int frame_id)
{
	auto size_in  = this->N;
	auto size_out = this->N_mod;
	auto bps      = this->bits_per_symbol;

	auto main_loop_size = size_in / bps;
	for (auto i = 0; i < main_loop_size; i++)
	{
		// determine the symbol with a lookup table
		unsigned idx = 0;
		for (auto j = 0; j < bps; j++)
			idx += unsigned(unsigned(1 << j) * X_N1[i * bps +j]);
		auto symbol = cstl->get_real(idx);

		X_N2[i] = symbol;
	}

	// last elements if "size_in" is not a multiple of the number of bits per symbol
	if (main_loop_size * bps < size_in)
	{
		unsigned idx = 0;
		for (auto j = 0; j < size_in - (main_loop_size * bps); j++)
			idx += unsigned(unsigned(1 << j) * X_N1[main_loop_size * bps +j]);
		auto symbol = cstl->get_real(idx);

		X_N2[size_out -1] = symbol;
	}
}

template <typename B,typename R, typename Q, tools::proto_max<Q> MAX>
void Modem_flash<B,R,Q,MAX>
::_demodulate(const Q *Y_N1, Q *Y_N2, const int frame_id)
{
	if (!std::is_same<R,Q>::value)
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, "Type 'R' and 'Q' have to be the same.");

	if (!std::is_floating_point<Q>::value)
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, "Type 'Q' has to be float or double.");

	if (!this->n->is_set())
		throw tools::runtime_error(__FILE__, __LINE__, __func__, "No noise has been set");

	auto size = this->N;
	auto thresholds_per_symbol = thresholder->get_n_thresholds_per_symbol();
	
	//Update thresholds for specific noise
	thresholder->update_thresholds();

	//Loop over noised symbols, check and interpret thresholds
	std::vector<R> readout(thresholds_per_symbol);
	for (auto i_symbol = 0; i_symbol < nbr_symbols; i_symbol++){
		for (auto i_threshold = 0; i_threshold < thresholds_per_symbol; i_threshold++){
			
			readout[i_threshold] = Y_N1[i_symbol] < thresholder->get_threshold(i_threshold) ? 0 : 1;
			Y_N2[i_symbol] = thresholder->interpret_readout(readout);
		}
	}
	
}


}
}
