#include <gtest/gtest.h>
#include <aff3ct_extension.hpp>

#include <string> 

namespace aff3ct
{
    namespace tools
    {

        TEST(FlashReaderTest, slc)
        {
            std::string fpath("test_data/reader");
            Flash_reader<float, float> slc(Flash_reader<float, float>::lower, Flash_reader<float, float>::hard, fpath);
        }

    }
}