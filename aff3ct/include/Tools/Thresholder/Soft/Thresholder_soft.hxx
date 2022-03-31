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
Thresholder_user<R>
::Thresholder_user(const std::string& const_path)
: Thresholder<R>(read_Thresholder(const_path), "User<C>")
{
}

template <typename R>
std::vector<typename Thresholder_user<R>::S> Thresholder_user<R>
::read_Thresholder(const std::string& const_path)
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

	std::vector<S> Thresholder;
	auto sqrt_es = (R)0;
	std::string temp;
	while (std::getline(const_file, temp))
	{
		if (temp[0] == '#') continue;

		std::istringstream buffer(temp);
		std::vector<R> line((std::istream_iterator<R>(buffer)), std::istream_iterator<R>());

		if (line.size() >= 3)
		{
			std::stringstream message;
			message << "'line.size()' has to be smaller than 3 ('line.size()' = " << line.size() << ").";
			throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
		}

		if (line.size() == 2)
			Thresholder.push_back(S(line[0],line[1]));
		else
			Thresholder.push_back(S(line[0]));

		sqrt_es += std::norm(Thresholder.back());
	}
	sqrt_es = std::sqrt(sqrt_es/Thresholder.size());

	for (unsigned i = 0; i < Thresholder.size(); i++)
		Thresholder[i] /= S(sqrt_es);

	return Thresholder;
}

}
}
