/*
 * \file
 * \brief Class tools::Flash_cell
 */
#ifndef FLASH_CELL_HPP__
#define FLASH_CELL_HPP__

#include <vector>
#include <aff3ct.hpp>

namespace aff3ct
{
namespace tools
{
template <typename R>
class Flash_cell
{
public:
	explicit Flash_cell(const int cell_type);

	R get_level(const unsigned int level_index);

	//Reads n_pages of data from Q and returns the mapped level index
	unsigned int get_level_index(const int symbol);

	unsigned int get_symbol(const unsigned int level_index);
	std::vector<unsigned int>& get_threshold_indexes(const unsigned int page_type);

	unsigned int get_n_pages();
	unsigned int get_n_levels();
	unsigned int get_n_threshold_indexes();

	enum page_type {lower, upper, extra};
	enum cell_type {SLC = 2, MLC = 4, TLC = 8};

private:
	
	unsigned int n_threshold_indexes;
	unsigned int n_levels;
	unsigned int n_pages;

	std::vector<std::vector<unsigned int>> threshold_indexes;

	std::vector<unsigned int> level_index_map;
	std::vector<unsigned int> symbol_map;
	std::vector<unsigned int> init_gray(const int cell_type);
};

}
}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "Tools/Flash/Cell/Flash_cell.hxx"
#endif

#endif // FLASH_CELL_HPP__