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
template <typename R>
Flash_cell<R>
::Flash_cell(const int cell_type) :
	n_threshold_indexes(cell_type - 1),
	n_levels(cell_type),
	n_pages(std::log2(cell_type)),

	level_index_map(cell_type),
	symbol_map(cell_type)
{
	
	//Could be automated
	switch (cell_type){
		case SLC:
			this -> threshold_indexes = {{0}};
		case MLC:
			this -> threshold_indexes = {{1},{0,2}};
		case TLC:
			this -> threshold_indexes = {{3}, {1,5}, {0, 2, 4, 6}};
		default:
			//Throw error
			break;
	}
}

template <typename R>
std::vector<unsigned int> Flash_cell<R>
::init_gray(const int cell_type){	
	int gray_i;
	for (int i = 0; i < cell_type; i++){
		gray_i = i^(i>>1);

		this->level_index_map[i] = gray_i;
		this->symbol_map[gray_i] = i;
	}
}

template <typename R>
unsigned int Flash_cell<R>
::get_level_index(const int symbol)
{
	return this->level_index_map[symbol];
}

template <typename R>
unsigned int  Flash_cell<R>
::get_symbol(const unsigned int level_index)
{
	return this->symbol_map[level_index];
}

template <typename R>
std::vector<unsigned int>&  Flash_cell<R>
::get_threshold_indexes(const unsigned int page_type)
{
	return this->threshold_indexes[page_type];
}

template <typename R>
unsigned int Flash_cell<R>
::get_n_pages(){
	return this->n_pages;
}

template <typename R>
unsigned int Flash_cell<R>
::get_n_levels()
{
	return this->level_index_map.size();
}

template <typename R>
unsigned int Flash_cell<R>
::get_n_threshold_indexes()
{
	return this->n_threshold_indexes;
}

}
}
