#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <fstream>
#include <cmath>
#include <iterator>
#include <sstream>

#include <aff3ct.hpp>

namespace aff3ct
{
namespace tools
{
Flash_cell
::Flash_cell(const unsigned cell_type) :
	n_threshold_indexes(cell_type - 1),
	n_levels(cell_type)
{
	n_pages = (unsigned)std::log2(cell_type);

	level_index_map = std::vector<unsigned>(n_levels);
	symbol_map = std::vector<unsigned>(n_levels);
	init_gray(cell_type);

	//Could be automated
	switch (cell_type){
		case SLC:
			threshold_indexes.push_back({0});
		case MLC:
			threshold_indexes.push_back({1});
			threshold_indexes.push_back({0,2});
		case TLC:
			threshold_indexes.push_back({3});
			threshold_indexes.push_back({1,5});
			threshold_indexes.push_back({0,2,4,6});
		default:
			//Throw error
			break;
	}
}


void Flash_cell
::init_gray(const unsigned cell_type){	
	unsigned gray_i;
	for (auto i = 0; i < (unsigned)cell_type; i++){
		gray_i = i^(i>>1);

		level_index_map[i] = gray_i;
		symbol_map[gray_i] = i;
	}
}

unsigned Flash_cell
::get_level_index(const unsigned symbol)
{
	return level_index_map[symbol];
}

unsigned  Flash_cell
::get_symbol(const unsigned level_index)
{
	return symbol_map[level_index];
}

std::vector<unsigned> Flash_cell
::get_threshold_indexes(const unsigned page_type)
{
	return threshold_indexes[page_type];
}

unsigned Flash_cell
::get_n_pages() const
{
	return n_pages;
}

unsigned Flash_cell
::get_n_levels() const
{
	return level_index_map.size();
}

unsigned Flash_cell
::get_n_threshold_indexes() const
{
	return n_threshold_indexes;
}

}
}
