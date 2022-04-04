/*! COPIED FROM CONSTELLATION FOR REFERENCE*/
#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <utility>
#include <sstream>
#include <cmath>
#include <aff3ct.hpp>

namespace aff3ct
{
namespace tools
{
template <typename R>
Thresholder<R>
::Thresholder(const unsigned n_tps, const std::string& name)
: n_tps(n_tps),
  name(name)
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
::get_n_thresholds_per_symbol() const
{
	return n_tps;
}

template <typename R>
const R Thresholder<R>
::operator[](const size_t idx) const
{
	return thresholds[idx];
}

template <typename R>
const R Thresholder<R>
::get_threshold(const size_t idx) const
{
	return thresholds[idx];
}

template <typename R>
R Thresholder<R>
::interpret_readout(const R* read_start, const R* read_stop) const
{
	throw tools::unimplemented_error(__FILE__, __LINE__, __func__);
}

template <typename R>
void Thresholder<R>
::update_thresholds() const
{
	throw tools::unimplemented_error(__FILE__, __LINE__, __func__);
}

}
}
