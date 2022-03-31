/*
 * \file
 * \brief Class tools::Thresholder
 */
#ifndef THRESHOLDER_HPP__
#define THRESHOLDER_HPP__

#include <cstdint>
#include <cstddef>
#include <string>
#include <aff3ct.hpp>

namespace aff3ct
{
namespace tools
{
/*
 * \brief describe a thresholder
 */
template <typename R>
class Thresholder
{
public:

	virtual ~Thresholder() = default;

	/*
	 * \brief get the name of the thresholder
	 */
	const std::string& get_name() const;

	/*
	 * \brief get the number of thresholds to test per symbol
	 */
	unsigned get_n_thresholds_per_symbol() const;

	/*
	 * \param idx is the index of the wanted threshold.
	 * \return the threshold level of the thresholder at the given index
	 */
	const R operator[] (const size_t idx) const;
	const R get_threshold (const size_t idx) const;

	/*
	 * \param read_start is a pointer to the start of readout
	 * \param read_stop is a pointer to the end of the readout
	 * \return a demodulated value interpreted from the readout
	 */ 
	virtual R interpret_readout(const Q *read_start, const Q *read_stop);

protected:
	/*
	 * \param n_tps is the number of thresholds per symbol
	 * \param name is the name of the thresholder
	 */
	Thresholder(const unsigned n_tps, const std::string& name);

private:
	const unsigned n_tps;  // the number of thresholds per symbol
	const std::string name; // the name of the thresholder
	std::vector<R> thresholder; // The tresholds of the thresholder
};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Thresholder/Thresholder.hxx"
#endif

#endif // THRESHOLDER_HPP__