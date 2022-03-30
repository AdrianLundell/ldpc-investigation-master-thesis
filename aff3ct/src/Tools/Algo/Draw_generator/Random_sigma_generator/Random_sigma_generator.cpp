#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct::tools;

template <typename R>
Random_sigma_generator<R>
::Random_sigma_generator(const int seed)
: Random_sigma_generator<R>()
{
	this->set_seed(seed);
}

template <typename R>
void Random_sigma_generator<R>
::set_seed(const int seed)
{
	rd_engine.seed(seed);
}

template <typename R>
void Random_sigma_generator<R>
::generate(std::map<R,R> &sigma_map, const unsigned size, const R sigma_tot)
{
	normal_dist = std::normal_distribution<R>(mu, sigma);

	for (unsigned i = 0; i < length; i++)
		noise[i] = this->normal_dist(this->rd_engine);
}

// ==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::tools::Random_sigma_generator<R_32>;
template class aff3ct::tools::Random_sigma_generator<R_64>;
#else
template class aff3ct::tools::Random_sigma_generator<R>;
#endif
// ==================================================================================== explicit template instantiation
