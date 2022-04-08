#include <gtest/gtest.h>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

namespace aff3ct 
{
namespace tools
{

// template <typename R>
// class FlashCellFixture : public ::testing::Test {
//     protected:
//         void SetUp() override {
//             slc_cell = new Flash_cell<R>(Flash_cell<R>::SLC);
//             mlc_cell = new Flash_cell<R>(Flash_cell<R>::MLC);
//             tlc_cell = new Flash_cell<R>(Flash_cell<R>::TLC);
//         }

//         Flash_cell<R> slc_cell;
//         Flash_cell<R> mlc_cell;
//         Flash_cell<R> tlc_cell;
// };

// //See types.h
// using rTypes = ::testing::Types<float, double>;
// TYPED_TEST_SUITE(FlashCellFixture, rTypes);

// Demonstrate some basic assertions.
TEST(FlashCellFixture, slc)
{
    Flash_cell<float> slc_cell(Flash_cell<float>::SLC);
    Flash_cell<float> mlc_cell(Flash_cell<float>::MLC);
    Flash_cell<float> tlc_cell(Flash_cell<float>::TLC);
}




}
}