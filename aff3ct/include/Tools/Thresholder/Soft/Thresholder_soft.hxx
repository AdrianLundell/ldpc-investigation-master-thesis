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
::update_thresholds(){
	float x = 5.f; //Should be calculated from noise
	float y = 5.f; //Should be calculated from noisetypedef vector <int> ints_t;

	float x1, x2, y1, y2;
	//Loop over data to find boundary points assuming increasing values row wise and (snr, ratio) within limits
	for (auto i = 0; i < this->data.size(); i++)
	{
		if ((this->data)[i][0] >= x && this->data[i][1] >= y)
		{
			//Weighted mean interpolation: https://en.wikipedia.org/wiki/Bilinear_interpolation#Weighted_mean
			std::vector<R> q22 = this->data[i];
			std::vector<R> q21 = this->data[i - 1];
			std::vector<R> q12 = this->data[i - this->n_snrs];
			std::vector<R> q11 = this->data[i - this->n_snrs - 1];

			x1 = q11[0];
			x2 = q22[0];
			y1 = q11[1];
			y2 = q22[1];
			
			float w11, w12, w21, w22, c;
			c = (x2 - x1)*(y2 - y1);
			w11 = (x2 - x)*(y2 - y)/c;
			w12 = (x2 - x)*(y - y1)/c;
			w21 = (x - x1)*(y2 - y)/c;
			w22 = (x - x1)*(y - y1)/c;

			int j;
			for (auto i = 0; i < this->n_thresholds; i++)
			{
				j = i + 2;
				thresholds[i] = w11*q11[j] + w21*q21[j] + w12*q12[j] * w22*q22[j];	
			}

			for (auto i = 0; i < this->n_llrs; i++)
			{
				j = i + 2 + this->n_thresholds;
				llrs[i] = w11*q11[j] + w21*q21[j] + w12*q12[j] * w22*q22[j];	
			}
			
			break;
		}
	}

	
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

	bool is_unique;
	for (auto x1 : x){
		is_unique = true;
		for (auto x2 : unique_elements)
			is_unique = (std::abs(x1 - x2) > 0.0001) && is_unique;
		
		if (is_unique)
			unique_elements.push_back(x1);
	}

	return unique_elements.size();
}

}
}
