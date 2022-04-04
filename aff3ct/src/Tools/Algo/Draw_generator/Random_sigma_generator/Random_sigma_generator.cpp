// #include <aff3ct.hpp>
// #include <aff3ct_extension.hpp>

// #include <set>

// using namespace aff3ct::tools;

// template <typename R>
// Random_sigma_generator<R>::Random_sigma_generator(const int seed)
// 	: Random_sigma_generator<R>()
// {
// 	this->set_seed(seed);
// }

// template <typename R>
// void Random_sigma_generator<R>::set_seed(const int seed)
// {
// 	rd_engine.seed(seed);
// }

// template <typename R>
// void Random_sigma_generator<R>::generate(std::map<R, R> &sigma_map, const R sigma_tot, const R min_sigma)
// {
// 	n_bins = (int)sigma_tot / min_sigma;
// 	uniform_int_dist = std::uniform_int_distribution<int>(1, n_bins - 1);

// 	std::set<int> rand_bin_levels;
// 	unsigned i = 0;
// 	unsigned max_i = sigma_map.size() - 1;
// 	while (i < max_i)
// 	{
// 		int rand_bin_level = this->uniform_int_dist(this->rd_engine);
// 		if (rand_bin_levels.find(rand_bin_level) == rand_bin_levels.end())
// 		{
// 			rand_bin_levels.insert(rand_bin_level);
// 			i++;
// 		}
// 	}
// 	rand_bin_levels.insert(n_bins);

// 	auto bin_levels_it = random_bin_levels.begin();
// 	int previous_bin = 0;
// 	for (auto sigma_it = sigma_map.begin(); sigma_it != sigma_it.end(); sigma_it++)
// 	{
// 		sigma_it->second = (R)(*bin_levels_it - previous_bin) * min_sigma;
// 		previous_bin = *bin_levels_it;
// 		bin_levels_it++;
// 	}
// }

// // ==================================================================================== explicit template instantiation
// #include "Tools/types.h"
// #ifdef AFF3CT_MULTI_PREC
// template class aff3ct::tools::Random_sigma_generator<R_32>;
// template class aff3ct::tools::Random_sigma_generator<R_64>;
// #else
// template class aff3ct::tools::Random_sigma_generator<R>;
// #endif
// // ==================================================================================== explicit template instantiation
