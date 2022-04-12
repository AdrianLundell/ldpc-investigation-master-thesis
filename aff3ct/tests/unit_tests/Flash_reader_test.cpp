#include <gtest/gtest.h>
#include <aff3ct_extension.hpp>

#include <string> 
#include <vector>

namespace aff3ct
{
namespace tools
{

TEST(FlashReaderTest, static_hard_slc)
{
    float sigma_tot = 1;
    float sigma_min = 0.5;
    unsigned n_sigmas = 2;
    Sigma_asymmetric<float> s;
    s.set_sigmas(sigma_tot, n_sigmas, sigma_min, 4, 4);
    s.generate_sigmas();

    std::vector<float> voltage_levels = {-1,1};
    aff3ct::module::Channel_AWGN_asymmetric<float> c(1, voltage_levels, s, 0, 1);
    c.set_noise(s);

    Flash_reader<float, float> slc(Flash_reader<float, float>::lower, Flash_reader<float, float>::hard, "test_data/reader_static_hard.txt");
    slc.update(c);

    //For symmetric sigmas
    //Extreme level differences
    ASSERT_FLOAT_EQ(slc.read(-10000, {0}), 0);
    ASSERT_FLOAT_EQ(slc.read(10000, {0}), 1);
    //Normal level differences
    ASSERT_FLOAT_EQ(slc.read(-1, {0}), 0);
    ASSERT_FLOAT_EQ(slc.read(1, {0}), 1);
    //Small voltage differences
    ASSERT_FLOAT_EQ(slc.read(-0.00001, {0}), 0);
    ASSERT_FLOAT_EQ(slc.read(+0.00001, {0}), 1);

    //Test same result regardless of sigmas
}

TEST(FlashReaderTest, static_soft_tlc)
{
    float sigma_tot = 1;
    float sigma_min = 0.5;
    unsigned n_sigmas = 2;
    Sigma_asymmetric<float> s;
    s.set_sigmas(sigma_tot, n_sigmas, sigma_min, 4, 4);
    s.generate_sigmas();

    std::vector<float> voltage_levels = {-1,1};
    aff3ct::module::Channel_AWGN_asymmetric<float> c(1, voltage_levels, s, 0, 1);
    c.set_noise(s);

    Flash_reader<float, float> slc(Flash_reader<float, float>::lower, Flash_reader<float, float>::soft_single, "test_data/reader_static_soft.txt");
    slc.update(c);

    //For symmetric sigmas
    ASSERT_FLOAT_EQ(slc.read(-1, {0}), 10);
    ASSERT_FLOAT_EQ(slc.read(-0.1, {0}), 20);
    ASSERT_FLOAT_EQ(slc.read(0.1, {0}), 30);
    ASSERT_FLOAT_EQ(slc.read(1, {0}), 40);
}


}
}