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
Thresholder_soft<R>
::Thresholder_soft(const std::string& const_path)
: Thresholder<R>(3, "soft<R>")
{
	thresholds = std::vector<R>(3);
	llrs = std::vector<R>(4);
	data = std::vector<std::vector<R>>();

	init_data(const_path, data);
}

template <typename R>
R Thresholder_soft<R>
::interpret_readout(const std::vector<R> &readout){
	
	//Return LLR left of the first positive readout bit
	for (int i = 0; i < readout.size(); i++){
		if (readout[i] == 1)
			return llrs[i];
	}
	return llrs.back();
}

template <typename R>
void Thresholder_soft<R>
::update_thresholds(const tools::Noise<R>& noise){


	// float sigma_db = noise.get_sigmas();
	// float sigma_ratio = noise.get_sigmas();

	// thresholds = data[];
	// llrs = data[];

	// std::vector<R>& x0, x1, x2, x3; 

	// //Find references to 
	// for (auto const &x : data)
	// {
	// 	if (x[0] >= sigma_db && x[1] >= sigma_ratio)
	// 	{
	// 		x1 = x;
	// 		break;
	// 	}
	// 	x0 = x1;
	// }
}

template <typename R>
void Thresholder_soft<R>
::init_data(const std::string& const_path,
			std::vector<std::vector<R>>& data)
{
	if (const_path.empty())
		throw tools::invalid_argument(__FILE__, __LINE__, __func__, "'const_path' should not be empty.");

	std::ifstream const_file(const_path);

	if (const_file.fail())
	{
		std::stringstream message;
		message << "Opening 'const_path' (= " << const_path << ") has failed.";
		throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
	}


	std::vector<R> snrs = {};
	std::vector<R> ratios = {};
	data.clear();

	std::string temp;
	while (std::getline(const_file, temp))
	{
		if (temp[0] == '#') continue;

		std::istringstream buffer(temp);
		std::vector<R> line((std::istream_iterator<R>(buffer)), std::istream_iterator<R>());

		if (!line.size()== 9) //We want to read sigma_db, sigma_ratio, n_tps thresholds and n_tps+1 llrs 
		{
			std::stringstream message;
			message << "'line.size()' has to be 9 ('line.size()' = " << line.size() << ").";
			throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
		} else {
			data.push_back(line);
			snrs.push_back(line[0]);
			ratios.push_back(line[1]);
		}
	}

	this->n_snrs = count_unique(snrs);
	this->n_ratios = count_unique(ratios);
}

template <typename R>
int Thresholder_soft<R>
::count_unique(const std::vector<R>& x)
{
	std::vector<R> unique_elements = {x[0]};

	for (auto x1 : x){
		for (auto x2 : unique_elements) {
			if (std::abs(x1 - x2) > 0.0001)
			{
				unique_elements.push_back(x1);
			}
		}
	}

	return unique_elements.size();
}

}
}
