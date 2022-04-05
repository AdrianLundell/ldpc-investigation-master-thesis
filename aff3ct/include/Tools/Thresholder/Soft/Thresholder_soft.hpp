/*
 * \file
 * \brief Class tools::Thresholder_soft.
 */
#ifndef THRESHOLDER_soft_HPP__
#define THRESHOLDER_soft_HPP__

#include <vector>

#include <aff3ct.hpp>

namespace aff3ct
{
namespace tools
{
template <typename R>
class Thresholder_soft : public Thresholder<R>
{
public:

	/*
	 * \param threshold_path is the path to soft thresholds
	 */
	explicit Thresholder_soft(const std::string& threshold_path);

	void update_thresholds();

	R interpret_readout(const std::vector<R> &readout);

private:
 
	/*
	 * \brief load data from txt file
	 */
	void init_data(const std::string& const_path, 
					std::vector<std::vector<R>>& data);

	/*
	 * \bried count number of unique numbers, floating precision compatible
	 */
	int count_unique(const std::vector<R>& x);

	std::vector<R> llrs;
	std::vector<R> thresholds;
	std::vector<std::vector<R>> data;

	int n_snrs;
	int n_ratios;
};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Thresholder/Soft/Thresholder_soft.hxx"
#endif

#endif // THRESHOLDER_soft_HPP__