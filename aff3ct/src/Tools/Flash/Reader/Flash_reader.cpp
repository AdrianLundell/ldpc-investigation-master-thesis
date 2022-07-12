// #ifndef _USE_MATH_DEFINES
// #define _USE_MATH_DEFINES
// #endif
#include <fstream>
#include <cmath>
#include <iterator>
#include <sstream>

#include <aff3ct_extension.hpp>

namespace aff3ct
{
	namespace tools
	{
		template <typename R, typename Q>
		Flash_reader<R, Q>::Flash_reader(const unsigned page_type, const unsigned read_type, std::string fpath) : thresholds(page_type),
																												  bin_values(page_type),
																												  my_page_type(page_type),
																												  my_read_type(read_type)
		{
			n_thresholds = read_type;
			n_bin_values = read_type + 1;
			init_data(fpath);
		}

		template <typename R, typename Q>
		void Flash_reader<R, Q>::update(c &channel)
		{

			float x, y;
			for (auto i = 0; i < this->get_page_type(); i++)
			{
				float sigma_ave = (float)channel.get_sigma_ave();
				this->thresholds[i] = std::vector<R>(this->get_n_thresholds());
				this->bin_values[i] = std::vector<R>(this->get_n_bin_values());

				this->_update(sigma_ave, this->thresholds[i], this->bin_values[i]);
			}
		}

		template <typename R, typename Q>
		void Flash_reader<R, Q>::_update(const float sigma_ave, std::vector<R> &thresholds, std::vector<Q> &bin_values)
		{
			// Important! When generating file from compute dmc, make sure to match the skewness with the one used in aff3ct

			const float sigma = 2 * sigma_ave; // sigma represents a total sigma

			// Bisection search to find the interval of sigma to interpolat from
			unsigned lower_idx = 0;
			unsigned upper_idx = this->data.size() - 1;

			// Firstly check that sigma_tot is within range
			float sigma_lower = data[lower_idx][1];
			float sigma_upper = data[upper_idx][1];
			if (sigma < sigma_lower || sigma > sigma_upper)
			{
				std::stringstream message;
				message << "sigma is out of range. sigma = " << sigma << ". Should be in range [" << sigma_lower << "," << sigma_upper << "].";
				throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
			}

			// Search for interpolation range thorugh bisections search
			unsigned mid_idx = (upper_idx + lower_idx) / 2;
			float sigma_mid = data[mid_idx][1];

			while (upper_idx - lower_idx != 1)
			{
				if (sigma > sigma_mid)
					lower_idx = mid_idx;
				else
					upper_idx = mid_idx;

				mid_idx = (upper_idx + lower_idx) / 2;
				sigma_mid = data[mid_idx][1];
			}

			// Interpolate values
			sigma_lower = data[lower_idx][1];
			sigma_upper = data[upper_idx][1];
			float weight = (sigma - sigma_lower) / (sigma_upper - sigma_lower);

			unsigned j;
			for (auto k = 0; k < this->n_thresholds; k++)
			{
				j = k + 5;
				thresholds[k] = data[lower_idx][j] + weight * (data[upper_idx][j] - data[lower_idx][j]);
			}

			for (auto k = 0; k < this->n_bin_values; k++)
			{
				j = k + 5 + this->n_thresholds;
				bin_values[k] = data[lower_idx][j] + weight * (data[upper_idx][j] - data[lower_idx][j]);
			}
		}

		template <typename R, typename Q>
		void Flash_reader<R, Q>::init_data(const std::string &fpath)
		{
			// Reads a file describing discretized functions z_i = f_i(x, y)
			// with columns x, y, z_1, z_2, ... z_n
			// And counts the dimension of the grid (x,y)

			if (fpath.empty())
				throw tools::invalid_argument(__FILE__, __LINE__, __func__, "'fpath' should not be empty.");

			std::ifstream file(fpath);

			if (file.fail())
			{
				std::stringstream message;
				message << "Opening 'fpath' (= " << fpath << ") has failed.";
				throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
			}

			std::vector<R> x_vec = {};
			std::vector<R> y_vec = {};
			this->data = std::vector<std::vector<R>>();

			std::string temp;
			while (std::getline(file, temp))
			{
				if (temp[0] == '#')
					continue;

				std::istringstream buffer(temp);
				std::vector<R> line((std::istream_iterator<R>(buffer)), std::istream_iterator<R>());

				if (!line.size() == 2 + this->n_thresholds + this->n_bin_values)
				{
					std::stringstream message;
					message << "'line.size()' has to match read type ('line.size()' = " << line.size() << ").";
					throw tools::runtime_error(__FILE__, __LINE__, __func__, message.str());
				}
				else
				{
					data.push_back(line);
					x_vec.push_back(line[0]);
					y_vec.push_back(line[1]);
				}
			}

			this->n_x = count_unique(x_vec);
			this->n_y = count_unique(y_vec);
		}

		template <typename R, typename Q>
		int Flash_reader<R, Q>::count_unique(const std::vector<R> &x)
		{
			std::vector<R> unique_elements = {x[0]};
			bool is_unique;

			for (auto x1 : x)
			{
				is_unique = true;
				for (auto x2 : unique_elements)
					is_unique = (std::abs(x1 - x2) > 0.0001) && is_unique;

				if (is_unique)
					unique_elements.push_back(x1);
			}

			return unique_elements.size();
		}

		template <typename R, typename Q>
		Q Flash_reader<R, Q>::read(const R level, const std::vector<unsigned> &threshold_indexes)
		{
			for (auto i : threshold_indexes)
			{
				for (auto j = 0; j < this->get_n_thresholds(); j++)
				{
					if (this->get_threshold(i, j) > level)
						return get_bin_value(i, j);
				}
			}

			return bin_values.back().back();
		}

		template <typename R, typename Q>
		R Flash_reader<R, Q>::get_threshold(const unsigned threshold_index, const unsigned soft_index)
		{
			return this->thresholds[threshold_index][soft_index];
		}

		template <typename R, typename Q>
		Q Flash_reader<R, Q>::get_bin_value(const unsigned threshold_index, const unsigned bin_index)
		{
			return this->bin_values[threshold_index][bin_index];
		}

		template <typename R, typename Q>
		unsigned Flash_reader<R, Q>::get_read_type()
		{
			return this->my_read_type;
		}

		template <typename R, typename Q>
		unsigned Flash_reader<R, Q>::get_page_type()
		{
			return this->my_page_type;
		}

		template <typename R, typename Q>
		unsigned Flash_reader<R, Q>::get_n_thresholds()
		{
			return this->n_thresholds;
		}

		template <typename R, typename Q>
		unsigned Flash_reader<R, Q>::get_n_bin_values()
		{
			return this->n_bin_values;
		}
	}
}

//==================================================================================== explicit template instantiation
#include "Tools/types.h"
#ifdef AFF3CT_MULTI_PREC
template class aff3ct::tools::Flash_reader<float, float>;
// template class aff3ct::tools::Flash_reader<R_32, R_32>;
// template class aff3ct::tools::Flash_reader<R_64, R_64>;
// template class aff3ct::tools::Flash_reader<R_64, Q_32>;
// TODO: More template instanciations
#else
template class aff3ct::tools::Flash_reader<R, Q>;
#endif
//==================================================================================== explicit template instantiation
