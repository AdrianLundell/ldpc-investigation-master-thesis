#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>
#include <random>
#include <cmath>

using namespace aff3ct::tools;

template <typename R>
Random_sigma_generator<R>::Random_sigma_generator(const int seed)
	: Draw_generator<R>()
{
	this->set_seed(seed);
}

template <typename R>
void Random_sigma_generator<R>::set_seed(const int seed)
{
	rd_engine.seed(seed);
}

template <typename R>
void Random_sigma_generator<R>::generate(std::vector<R> &sigmas, const R sigma_tot, const R sigma_min)
{
	this->uniform_real_dist = std::uniform_real_distribution<R>((R)0.0, (R)1.0);
	unsigned n_sigmas = sigmas.size();
	std::vector<R> sigma_vars(n_sigmas, (R)0.0);
	R sigma_vars_norm = 0.0;
	for (R sigma_var : sigma_vars)
	{
		sigma_var = this->uniform_real_dist(this->rd_engine);
		sigma_vars_norm += sigma_var * sigma_var;
	}
	for (unsigned i = 0; i < n_sigmas; i++)
	{
		sigmas[i] = sigma_min + sigma_vars[i] / sigma_vars_norm * (sigma_tot - sqrt(n_sigmas) * sigma_min);
	}
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
