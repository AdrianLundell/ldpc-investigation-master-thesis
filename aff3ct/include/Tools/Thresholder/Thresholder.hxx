/*! COPIED FROM CONSTELLATION FOR REFERENCE*/
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <utility>
#include <sstream>
#include <cmath>

#include <aff3ct>

namespace aff3ct
{
namespace tools
{
template <typename R>
bool has_complex_symbols(const Thresholder<R>& cstl)
{
	for (unsigned i = 0; i < cstl.get_n_symbols(); i++)
	{
		if (cstl[i].imag() != (R)0.)
			return true;
	}

	return false;
}

template <typename R>
Thresholder<R>
::Thresholder(std::vector<S>&& symbols, const std::string& name)
: n_bps((unsigned int)std::log2(symbols.size())),
  n_symbs(1 << n_bps),
  name(std::to_string(n_symbs) + name),
  Thresholder(std::move(symbols)),
  is_cplx(true)
{
	if (Thresholder.size() != this->get_n_symbols())
	{
		std::stringstream message;
		message << "'symbols.size()' has to be a power of 2 ('symbols.size()' = "
		        << symbols.size() << ").";
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, message.str());
	}
}

template <typename R>
Thresholder<R>
::Thresholder(const unsigned n_bps, const std::string& name)
: Thresholder(std::move(std::vector<S>((size_t)((int64_t)1 << n_bps))), name)
{
}

template <typename R>
Thresholder<R>
::Thresholder(const std::vector<S>& symbols, const std::string& name)
: Thresholder(std::move(std::vector<S>(symbols)), name)
{
}

template <typename R>
const std::string& Thresholder<R>
::get_name() const
{
	return name;
}

template <typename R>
unsigned Thresholder<R>
::get_n_bits_per_symbol() const
{
	return n_bps;
}

template <typename R>
unsigned Thresholder<R>
::get_n_symbols() const
{
	return n_symbs;
}

template <typename R>
bool Thresholder<R>
::is_complex() const
{
	return is_cplx;
}

template <typename R>
const typename Thresholder<R>::S& Thresholder<R>
::operator[](const size_t idx) const
{
	return Thresholder[idx];
}

template <typename R>
const typename Thresholder<R>::S& Thresholder<R>
::get_symbol(const size_t idx) const
{
	return Thresholder[idx];
}

template <typename R>
R Thresholder<R>
::get_imag(const size_t idx) const
{
	return Thresholder[idx].imag();
}

template <typename R>
R Thresholder<R>
::get_real(const size_t idx) const
{
	return Thresholder[idx].real();
}

template <typename R>
typename Thresholder<R>::S Thresholder<R>
::bits_to_symbol(const uint8_t bits[]) const
{
	throw tools::unimplemented_error(__FILE__, __LINE__, __func__);
}

template <typename R>
void Thresholder<R>
::build()
{
	try
	{
		std::vector<uint8_t> bits(this->get_n_bits_per_symbol());

		for (unsigned j = 0; j < this->get_n_symbols(); j++)
		{
			for (unsigned l = 0; l < this->get_n_bits_per_symbol(); l++)
				bits[l] = (uint8_t)((j >> l) & 1);
			this->Thresholder[j] = this->bits_to_symbol(bits.data());
		}
	}
	catch(tools::unimplemented_error&)
	{} // Thresholder has been filled by another way

	is_cplx = has_complex_symbols(*this);
}

}
}
