/*
 * \file
 * \brief Class tools::Thresholder_static.
 */
#ifndef THRESHOLDER_STATIC_HPP__
#define THRESHOLDER_STATIC_HPP__

#include <vector>

#include <aff3ct.hpp>

namespace aff3ct
{
namespace tools
{
template <typename R>
class Thresholder_static : public Thresholder<R>
{
public:

	/*
	 * \param threshold_path is the path to static thresholds
	 */
	explicit Thresholder_static(const std::string& threshold_path);

	R interpret_readout(const Q *read_start, const Q *read_stop);

private:
	static std::vector<S> read_thresholds(const std::string& const_path);
};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Thresholder/Static/Thresholder_static.hxx"
#endif

#endif // THRESHOLDER_STATIC_HPP__