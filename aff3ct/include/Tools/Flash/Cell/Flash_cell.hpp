/*!
 * \file
 * \brief Class tools::Flash_cell
 */
#ifndef FLASH_CELL_HPP_
#define FLASH_CELL_HPP_

#include <vector>

namespace aff3ct
{
namespace tools
{
class Flash_cell
{
public:
	explicit Flash_cell(const unsigned cell_type);

	//Reads n_pages of data from Q and returns the mapped level index
	unsigned get_level_index(const unsigned symbol);

	unsigned get_symbol(const unsigned level_index);
	std::vector<unsigned> get_threshold_indexes(const unsigned page_type);

	unsigned get_n_pages() const;
	unsigned get_n_levels() const;
	unsigned get_n_threshold_indexes() const;

	enum page_type {lower, upper, extra};
	enum cell_type {SLC = 2, MLC = 4, TLC = 8};

private:
	
	unsigned n_threshold_indexes;
	unsigned n_levels;
	unsigned n_pages;

	std::vector<std::vector<unsigned>> threshold_indexes;

	std::vector<unsigned> level_index_map;
	std::vector<unsigned> symbol_map;
	void init_gray(const unsigned cell_type);
};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Flash/Cell/Flash_cell.hxx"
#endif

#endif // FLASH_CELL_HPP_