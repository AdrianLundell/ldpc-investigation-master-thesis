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
#include <stdlib.h>  

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

namespace aff3ct
{
namespace module
{
template <typename B, typename R, typename Q>
Modem_flash_page<B,R,Q>
::Modem_flash_page(const int N,
				   tools::Flash_cell &cell,
				   tools::Flash_reader<R,Q> &reader,
				   const tools::Noise<R>& noise,
				   const int n_frames)
: Modem<B,R,Q>(N,
               N, // N_mod
              noise,
              n_frames),
  bits_per_symbol(cell.get_n_pages()),
  nbr_symbols    (cell.get_n_levels()),
  cell(cell),
  reader(reader)
{
	const std::string name = "Modem_flash_page"; //TODO: Use naming info from cell/ reader
	this->set_name(name);

	// if (cell == nullptr)
		// throw tools::invalid_argument(__FILE__, __LINE__, __func__, "No cell given ('cell' = nullptr).");


	// if (reader == nullptr)
		// throw tools::invalid_argument(__FILE__, __LINE__, __func__, "No reader given ('reader' = nullptr).");
}

template <typename B, typename R, typename Q>
void Modem_flash_page<B,R,Q>
::set_noise(module::Channel_AWGN_asymmetric<R> &channel)
{
	//Kept from old code, unclear why these methdos are needed.
	//Modem<B,R,Q>::set_noise(channel.current_noise());
	//this->n->is_of_type_throw(tools::Noise_type::SIGMA);

	//Update the thresholds and bin values of the reader
	reader.update(channel);
}

template <typename B,typename R, typename Q>
void Modem_flash_page<B,R,Q>
::_modulate(const B *X_N1, R *X_N2, const int frame_id)
{

	for (auto i = 0; i < this->N; i++)
	{
		// determine the symbol with a lookup table by converting binary to decimal representation
		// generate random data not relevant for the page to be read
		unsigned idx = 0;
		for (auto j = 0; j < bits_per_symbol; j++)
			if (j == std::log2(reader.get_page_type()))
			{	
				idx += unsigned(unsigned(1 << j) * X_N1[i]);
			} else {
				idx += unsigned(unsigned(1 << j) * std::rand()%2);
			} 
		auto symbol = this->cell.get_level_index(idx);

		X_N2[i] = symbol;
	}
}

template <typename B,typename R, typename Q>
void Modem_flash_page<B,R,Q>
::_demodulate(const Q *Y_N1, Q *Y_N2, const int frame_id)
{
	if (!std::is_same<R,Q>::value)
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, "Type 'R' and 'Q' have to be the same.");

	if (!std::is_floating_point<Q>::value)
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, "Type 'Q' has to be float or double.");

	if (!this->n->is_set())
		throw tools::runtime_error(__FILE__, __LINE__, __func__, "No noise has been set");
	
	for (auto i = 0; i < this->N; i++){		
		Y_N2[i] = reader.read(Y_N1[i], cell.get_threshold_indexes(std::log2(reader.get_page_type())));
	}
}

}
}
