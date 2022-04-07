/*
 * \file
 * \brief Class tools::Flash_reader
 */
#ifndef FLASH_READER_HPP__
#define FLASH_READER_HPP__

#include <vector>

#include <aff3ct.hpp>

namespace aff3ct
{
namespace tools
{
template <typename R, typename Q>
class Flash_reader : Flash_reader<R, Q>
{
public:

	explicit Flash_reader(const std::string& fpath);

	void update(const tools::Noise);

	Q read(const R level, const std::vector<unsigned>& threshold_indexes);

	enum read_type {hard = 1, soft_single = 3, soft_double = 5}

private:

	void init_data(const std::string& fpath);
	int count_unique(const std::vector<R>& x);

	float calculate_snr();
	float calculate_ratio();

	unsigned get_n_thresholds();
	R get_threshold(const unsigned threshold_index, const unsigned soft_index);
	Q get_bin_value(const unsigned threshold_index, const unsigned bin_index);

	std::vector<std::vector<Q>> bin_values;
	std::vector<std::vector<R>> thresholds;
	std::vector<std::vector<R>> data;

	int n_x;
	int n_y;
};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Thresholder/Soft/Thresholder_soft.hxx"
#endif

#endif // THRESHOLDER_soft_HPP__