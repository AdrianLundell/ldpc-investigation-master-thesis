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
		Flash_reader<R, Q>::Flash_reader(const int page_type, const int read_type, std::string fpath) : n_thresholds(read_type - 1),
																											  n_bin_values(read_type),
																											  thresholds(page_type),
																											  bin_values(page_type),
																											  my_page_type(page_type),
																											  my_read_type(read_type)
		{
			init_data(fpath);
		}

		template <typename R, typename Q>
		void Flash_reader<R, Q>::update(c &channel)
		{

			float x, y;
			for (auto i = 0; i < this->get_page_type(); i++)
			{
				x = channel.get_snr(i);
				y = channel.get_sigma_ratio(i);
				this->thresholds[i] = std::vector<R>(this->get_n_thresholds());
				this->bin_values[i] = std::vector<R>(this->get_n_bin_values());

				this->_update(x, y, this->thresholds[i], this->bin_values[i]);
			}
		}

		template <typename R, typename Q>
		void Flash_reader<R, Q>::_update(const float x, const float y, std::vector<R> &thresholds, std::vector<Q> &bin_values)
		{

			float x1, x2, y1, y2;
			// Loop over data to find boundary points assuming increasing values row wise and (snr, ratio) within limits
			for (auto i = 0; i < this->data.size(); i++)
			{
				if ((this->data)[i][0] >= x && this->data[i][1] >= y)
				{
					// Weighted mean interpolation: https://en.wikipedia.org/wiki/Bilinear_interpolation#Weighted_mean
					std::vector<R> &q22 = this->data[i];
					std::vector<R> &q21 = this->data[i - 1];
					std::vector<R> &q12 = this->data[i - this->n_x];
					std::vector<R> &q11 = this->data[i - this->n_x - 1];

					x1 = q11[0];
					x2 = q22[0];
					y1 = q11[1];
					y2 = q22[1];

					float w11, w12, w21, w22, c;
					c = (x2 - x1) * (y2 - y1);
					w11 = (x2 - x) * (y2 - y) / c;
					w12 = (x2 - x) * (y - y1) / c;
					w21 = (x - x1) * (y2 - y) / c;
					w22 = (x - x1) * (y - y1) / c;

					int j;
					for (auto i = 0; i < this->n_thresholds; i++)
					{
						j = i + 2;
						thresholds[i] = w11 * q11[j] + w21 * q21[j] + w12 * q12[j] + w22 * q22[j];
					}

					for (auto i = 0; i < this->n_bin_values; i++)
					{
						j = i + 2 + this->n_thresholds;
						bin_values[i] = w11 * q11[j] + w21 * q21[j] + w12 * q12[j] + w22 * q22[j];
					}

					break;
				}
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
					if (this->get_threshold(i, j) < level)
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
		int Flash_reader<R, Q>::get_read_type()
		{
			return this->my_read_type;
		}

		template <typename R, typename Q>
		int Flash_reader<R, Q>::get_page_type()
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
//template class aff3ct::tools::Flash_reader<R_32, R_32>;
//template class aff3ct::tools::Flash_reader<R_64, R_64>;
//template class aff3ct::tools::Flash_reader<R_64, Q_32>;
//TODO: More template instanciations
#else
template class aff3ct::tools::Flash_readerc<R, Q>;
#endif
//==================================================================================== explicit template instantiation
