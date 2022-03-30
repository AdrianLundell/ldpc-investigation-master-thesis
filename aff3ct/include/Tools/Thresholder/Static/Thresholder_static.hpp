/*!
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
	using typename Thresholder<R>::S;

	/*
	 * \param n_bps is the number of bits per symbol
	 */
	explicit Thresholder_static(const std::string& const_path);

private:
	static std::vector<S> read_thresholder(const std::string& const_path);
};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Thresholder/Static/Thresholder_static.hxx"
#endif

#endif // THRESHOLDER_STATIC_HPP__