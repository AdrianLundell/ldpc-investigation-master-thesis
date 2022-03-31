/*
 * \file
 * \brief Class tools::Constellation_user.
 */
#ifndef CONSTELLATION_FLASH_HPP__
#define CONSTELLATION_FLASH_HPP__

#include <vector>

#include <aff3ct.hpp>

namespace aff3ct
{
namespace tools
{
template <typename R>
class Constellation_flash : public Constellation<R>
{
public:
	using typename Constellation<R>::S;

	/*
	 * \param n_bps is the number of bits per symbol
	 */
	explicit Constellation_flash(const std::string& const_path);

private:
	static std::vector<S> read_constellation(const std::string& const_path);
};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Constellation/Flash/Constellation_flash.hxx"
#endif

#endif // CONSTELLATION_FLASH_HPP__