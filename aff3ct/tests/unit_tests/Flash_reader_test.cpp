#include <gtest/gtest.h>
#include <aff3ct_extension.hpp>

#include <cmath>
#include <string> 
#include <vector>

namespace aff3ct
{
namespace tools
{

TEST(FlashReaderTest, StaticHard)
{
    Flash_reader<float, float> reader(Flash_reader<float, float>::lower, Flash_reader<float, float>::hard, "test_data/reader_static_hard.txt");

    //Set symmetric noise
    Sigma_asymmetric<float> s1;
    std::vector<R> sigmas1 = {1, 1};
    s1.set_sigmas(sigmas1);
    std::vector<R> levels1 = {-1, 1};
    aff3ct::module::Channel_AWGN_asymmetric<float> c1(1, levels1, s1, 0, 1);
    c1.set_noise(s1);
    reader.update(c1);

    ASSERT_FLOAT_EQ(reader.read(-10000, {0}), 0);
    ASSERT_FLOAT_EQ(reader.read(10000, {0}), 1);
    ASSERT_FLOAT_EQ(reader.read(-1, {0}), 0);
    ASSERT_FLOAT_EQ(reader.read(1, {0}), 1);
    ASSERT_FLOAT_EQ(reader.read(-0.00001, {0}), 0);
    ASSERT_FLOAT_EQ(reader.read(+0.00001, {0}), 1);

    //Set asymmetric noise
    Sigma_asymmetric<float> s2;
    std::vector<R> sigmas2 = {1, 1};
    s2.set_sigmas(sigmas2);
    std::vector<R> levels2 = {-1, 1};
    aff3ct::module::Channel_AWGN_asymmetric<float> c2(1, levels2, s2, 0, 1);
    c2.set_noise(s2);
    reader.update(c2);

    ASSERT_FLOAT_EQ(reader.read(-10000, {0}), 0);
    ASSERT_FLOAT_EQ(reader.read(10000, {0}), 1);
    ASSERT_FLOAT_EQ(reader.read(-1, {0}), 0);
    ASSERT_FLOAT_EQ(reader.read(1, {0}), 1);
    ASSERT_FLOAT_EQ(reader.read(-0.00001, {0}), 0);
    ASSERT_FLOAT_EQ(reader.read(+0.00001, {0}), 1);
}

TEST(FlashReaderTest, DynamicHard)
{
    Flash_reader<float, float> reader(Flash_reader<float, float>::lower, Flash_reader<float, float>::hard, "test_data/reader_dynamic_hard.txt");

    //Set symmetric noise
    Sigma_asymmetric<float> s1;
    std::vector<R> sigmas1 = {1, 1};
    s1.set_sigmas(sigmas1);
    std::vector<R> levels1 = {-1, 1};
    aff3ct::module::Channel_AWGN_asymmetric<float> c1(1, levels1, s1, 0, 1);
    c1.set_noise(s1);
    reader.update(c1);

    //Expected interpolated threshold
    auto threshold = 50;

    ASSERT_FLOAT_EQ(reader.read(threshold + 0.01, {0}), 0);
    ASSERT_FLOAT_EQ(reader.read(threshold - 0.01, {0}), 1);
}

TEST(FlashReaderTest, static_soft_tlc)
{
    Sigma_asymmetric<float> s;
    std::vector<float> sigmas = {1,1};
    s.set_sigmas(sigmas);

    std::vector<float> voltage_levels = {-1,1};
    aff3ct::module::Channel_AWGN_asymmetric<float> c(1, voltage_levels, s, 0, 1);
    c.set_noise(s);

    Flash_reader<float, float> reader(Flash_reader<float, float>::lower, Flash_reader<float, float>::soft_single, "test_data/reader_static_soft.txt");
    reader.update(c);

    //For symmetric sigmas
    ASSERT_FLOAT_EQ(reader.read(-1, {0}), 10);
    ASSERT_FLOAT_EQ(reader.read(-0.1, {0}), 20);
    ASSERT_FLOAT_EQ(reader.read(0.1, {0}), 30);
    ASSERT_FLOAT_EQ(reader.read(1, {0}), 40);
}


}
}