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
	data = std::map<float, std::vector<R>>();

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
::update_thresholds(float sigma_ratio){

	std::vector<R> y0;
	std::vector<R> y1;
	float x0; 
	float x1;

	std::vector<R> y(8);
	float x = sigma_ratio;

	//Find stored data points closest to sigma_ratio
	//TODO: Corner case handling
	for (auto it = data.begin(); it != data.end(); it++)
	{
		if (it->first > x)
		{
			x0 = it->first;
			y0 = it->second;
			it++;
			x1 = it->first;
			y1 = it->second;
			break;
		}
	}

	//Interpolate all values
	for (auto i = 0; i < y.size(); i++) {
    	y[i] = y0[i] + (x-x0)*(y1[i] - y0[i])/(x1-x0);
	}

	thresholds = std::vector<R>(y.begin(), y.begin() + 3);
	llrs = std::vector<R>(y.begin() + 3, y.end());
}

template <typename R>
void Thresholder_soft<R>
::init_data(const std::string& const_path,
			std::map<float, std::vector<R>> data)
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

	std::string temp;
	while (std::getline(const_file, temp))
	{
		if (temp[0] == '#') continue;

		std::istringstream buffer(temp);
		std::vector<R> line((std::istream_iterator<R>(buffer)), std::istream_iterator<R>());

		float sigma_ratio;
		if (!line.size()== 8) //We want to read one sigma_ratio, n_tps thresholds and n_tps+1 llrs 
		{
			std::stringstream message;
			message << "'line.size()' has to be 8 ('line.size()' = " << line.size() << ").";
			throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
		} else {
			//Store each line in a vector mapped to its sigma_ratio
			for (auto i = 1; i < line.size(); i++){
				sigma_ratio = line[0];
				//data[sigma_ratio][i-1] = line[i];
			}
		}
	}
}

}
}
