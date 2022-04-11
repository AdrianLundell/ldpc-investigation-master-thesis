#include <gtest/gtest.h>
#include <aff3ct_extension.hpp>

namespace aff3ct 
{
namespace tools
{

TEST(FlashReaderTest, slc)
{
    Flash_reader<float, float> slc(Flash_reader<float,float>::lower, Flash_reader<float, float>::hard, "test_data/reader");
}


}
}