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
: Thresholder<R>(read_Thresholder(const_path), "soft<R>")
{
	std::vector<R> thresholds(3);
	std::vector<R> llrs(4);
	std::vector<std::unordered_map<float, R>> threshold_data(3);
	std::vector<std::unordered_map<float, R>> llr_data(4);
	
	init_data(const_path, threshold_data, llr_data);
}

template <typename R>
R Thresholder_soft<R>
::interpret_readout(const Q *read_start, const Q *read_stop){
	
	for (auto i = 0; i < read_stop; i++){
		if (Q + i == 1)
			return llr_data[i]
	}
	return llr_data[3]
}

template <typename R>
void Thresholder_soft<R>
::update_thresholds(float sigma_ratio){

	for 

	data = datapoints<
}

template <typename R>
R thresholder_soft<R>
::interpolate

template <typename R>
void Thresholder_soft<R>
::init_data(const std::string& const_path,
			std::vector<std::unorderd_map<float,R>> threshold_data, 
			std::vector<std::unorderd_map<float,R>> llr_data)
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

		if (!line.size()== 8) //We want to read one sigma_ratio, n_tps thresholds and n_tps+1 llrs 
		{
			std::stringstream message;
			message << "'line.size()' has to be 3 ('line.size()' = " << line.size() << ").";
			throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
		}

		else
			sigma_ratio = R(line[0]);

			threshold_data[0].insert({sigma_ratio, line[1]});
			threshold_data[1].insert({sigma_ratio, line[2]});
			threshold_data[2].insert({sigma_ratio, line[3]});

			threshold_data[0].insert({sigma_ratio, line[4]});
			threshold_data[1].insert({sigma_ratio, line[5]});
			threshold_data[2].insert({sigma_ratio, line[6]});
			threshold_data[3].insert({sigma_ratio, line[7]});
	}
}

}
}
