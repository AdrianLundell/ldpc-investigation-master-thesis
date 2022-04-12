#include <gtest/gtest.h>
#include <aff3ct_extension.hpp>
#include <math.h>

#include <vector>

using namespace aff3ct;
class RandomSigmaGeneratorTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
    }

    tools::Random_sigma_generator<float> sigma_generator;
    // void TearDown() override {}
};

TEST_F(RandomSigmaGeneratorTest, generate)
{
    // Expect two strings not to be equal.
    float sigma_tot = 3.0;
    float sigma_min = 0.5;
    unsigned n_sigmas = 4;
    std::vector<float> sigmas(n_sigmas, 0.0f);

    // Test different seeds
    for (int seed = 0; seed < 20; seed++)
    {
        sigma_generator.set_seed(seed);
        sigma_generator.generate(sigmas, sigma_tot, sigma_min);

        for (float sigma : sigmas)
        {
            EXPECT_GE(sigma, sigma_min) << "Each individual sigma must greater or equal to sigma_min";
            float sigma_upper_bound = sigma_tot + sigma_min * (1 - sqrt(n_sigmas));
            EXPECT_LE(sigma, sigma_upper_bound) << "Each individual sigma must be smaller than (sigma_tot+sigma_min*(1-sqrt(n_sigmas))";
        }
    }
}
