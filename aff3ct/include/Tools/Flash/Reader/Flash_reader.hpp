/*
 * \file
 * \brief Class tools::Flash_reader
 */
#ifndef FLASH_READER_HPP_
#define FLASH_READER_HPP_

#include <vector>
#include <aff3ct.hpp>

namespace aff3ct
{
	namespace tools
	{
		template <typename R, typename Q>
		class Flash_reader
		{
		public:
			explicit Flash_reader(const int page_type, const int read_type, const std::string fpath);

			typedef aff3ct::module::Channel_AWGN_asymmetric<R> c;
			void update(c &channel);
			int get_read_type();
			int get_page_type();

			Q read(const R level, const std::vector<unsigned> &threshold_indexes);

			enum read_type
			{
				hard = 1,
				soft_single = 3,
				soft_double = 5
			};
			enum page_type
			{
				lower = 1,
				upper = 2,
				extra = 4
			};

		private:
			void _update(const float x, const float y, std::vector<R> &thresholds, std::vector<Q> &bin_values);

			void init_data(const std::string &fpath);
			int count_unique(const std::vector<R> &x);

			unsigned get_n_thresholds();
			unsigned get_n_bin_values();
			R get_threshold(const unsigned threshold_index, const unsigned soft_index);
			Q get_bin_value(const unsigned threshold_index, const unsigned bin_index);

			std::vector<std::vector<Q>> bin_values;
			std::vector<std::vector<R>> thresholds;
			std::vector<std::vector<R>> data;

			unsigned n_thresholds;
			unsigned n_bin_values;
			int n_x;
			int n_y;
			int my_page_type; // Change enumeration index
			int my_read_type;
		};

	}
}

// #ifndef DOXYGEN_SHOULD_SKIP_THIS
// #include "Tools/Flash/Reader/Flash_reader.hxx"
// #endif

#endif // FLASH_READER_HPP_