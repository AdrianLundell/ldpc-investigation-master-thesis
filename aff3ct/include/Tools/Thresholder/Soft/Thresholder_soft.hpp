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

	R interpret_readout(const Q *read_start, const Q *read_stop);

private:
 
	void init_data(const std::string& const_path, 
					std::vector<std::unordered_map<float, R>> threshold_data,
					std::vector<std::unordered_map<float, R>> llr_data);
	
	std::vector<S> llrs;
	std::vector<S> thresholds;
	std::vector<std::unorderd_map<float,R>> threshold_data;
	std::vector<std::unorderd_map<float,R>> llr_data;

};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Thresholder/Soft/Thresholder_soft.hxx"
#endif

#endif // THRESHOLDER_soft_HPP__