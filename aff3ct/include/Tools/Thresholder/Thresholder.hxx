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
::Thresholder(const unsigned n_bps, const std::string& name)
: Thresholder(std::move(std::vector<S>((size_t)((int64_t)1 << n_bps))), name)
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
const typename R& Thresholder<R>
::operator[](const size_t idx) const
{
	return thresholder[idx];
}

template <typename R>
const typename R& Thresholder<R>
::get_threshold(const size_t idx) const
{
	return thresholder[idx];
}

template <typename R>
typename R Thresholder<R>
::interpret_readout(const Q *read_start, const Q *read_stop) const
{
	throw tools::unimplemented_error(__FILE__, __LINE__, __func__);
}

}
}
