#include <gtest/gtest.h>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

namespace aff3ct 
{
namespace tools
{

TEST(FlashReaderTest, slc)
{
    Flash_reader slc = Flash_reader(Flash_reader::lower, Flash_reader::hard);
}


}
}