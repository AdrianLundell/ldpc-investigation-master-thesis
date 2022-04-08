#include <gtest/gtest.h>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

namespace aff3ct 
{
namespace tools
{

TEST(FlashCellTest, slc)
{
    Flash_cell slc = Flash_cell(Flash_cell::SLC);

    ASSERT_EQ(slc.get_level_index(0), 0);
    ASSERT_EQ(slc.get_level_index(1), 1);

}

TEST(FlashCellTest, mlc)
{
    Flash_cell mlc = Flash_cell(Flash_cell::MLC);

    ASSERT_EQ(mlc.get_level_index(0), 0);
    ASSERT_EQ(mlc.get_level_index(1), 1);
    ASSERT_EQ(mlc.get_level_index(2), 3);
    ASSERT_EQ(mlc.get_level_index(3), 2);
}



TEST(FlashCellTest, tlc)
{
    Flash_cell tlc(Flash_cell::TLC);

    ASSERT_EQ(tlc.get_level_index(0), 0);
    ASSERT_EQ(tlc.get_level_index(1), 1);
    ASSERT_EQ(tlc.get_level_index(2), 3);
    ASSERT_EQ(tlc.get_level_index(3), 2);
    ASSERT_EQ(tlc.get_level_index(4), 6);
    ASSERT_EQ(tlc.get_level_index(5), 7);
    ASSERT_EQ(tlc.get_level_index(6), 5);
    ASSERT_EQ(tlc.get_level_index(7), 4);


}

}
}