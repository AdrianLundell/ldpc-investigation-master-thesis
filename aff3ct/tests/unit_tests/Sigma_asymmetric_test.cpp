#include <gtest/gtest.h>
#include <aff3ct_extension.hpp>
#include <aff3ct.hpp>

#include <math.h>
#include <vector>

using namespace aff3ct;

class SigmaAsymmetricTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        ebn0 = 0.01f;
        R = 1;
        esn0 = tools::ebn0_to_esn0(ebn0, R);
        sigma_tot = tools::esn0_to_sigma(esn0);
        sigma_min = 0.05;
        n_sigmas = 8;
    }

    tools::Sigma_asymmetric<float> sigmas;
    float ebn0;
    float R;
    float esn0;
    float sigma_tot;
    float sigma_min;
    unsigned n_sigmas;
};

TEST_F(SigmaAsymmetricTest, initialize)
{
    EXPECT_FALSE(sigmas.has_sigmas()) << "sigmas.has_sigmas() should be false at intitialization.";
    EXPECT_ANY_THROW(sigmas.check()) << "sigmas.check() should throw exception when sigmas.set_sigmas() has not been called.";
    EXPECT_ANY_THROW(sigmas.generate_sigmas()) << "sigmas.generate_sigmas() should throw exception when sigmas.set_sigmas() has not been called";

    sigmas.set_sigmas(sigma_tot, n_sigmas, sigma_min, ebn0, esn0);
    EXPECT_TRUE(sigmas.has_sigmas()) << "sigmas.has_sigmas() should be true after sigmas.set_sigmas() has been called";
    EXPECT_NO_THROW(sigmas.check()) << "sigmas.check() should not throw exception when sigmas.set_sigmas() has been called.";
    EXPECT_NO_THROW(sigmas.generate_sigmas()) << "sigmas.generate_sigmas() should not throw exception when sigmas.set_sigmas() has been called";

    EXPECT_EQ(sigmas.get_sigmas().size(), n_sigmas) << "sigmas.get_sigmas() should return a vector of size n_sigmas.";
    EXPECT_EQ(sigmas.get_noise(), sigma_tot) << "sigmas.get_noise() should return sigma_tot.";

    // Make sure it throws exceptions if unallowed values for sigma_min and sigma_tot are set
    EXPECT_ANY_THROW(sigmas.set_sigmas(sigma_tot, n_sigmas, sigma_tot, ebn0, esn0)) << "sigma_min=sigma_tot should not be allowed.";
    EXPECT_ANY_THROW(sigmas.set_sigmas(0.0f, n_sigmas, sigma_min, ebn0, esn0)) << "sigma_tot=0 should not be allowed.";
    EXPECT_ANY_THROW(sigmas.set_sigmas(sigma_tot, n_sigmas, 0.0f, ebn0, esn0)) << "sigma_min=0 should not be allowed.";
}

TEST_F(SigmaAsymmetricTest, generateSigmas)
{
    // Only need to test that sigmas.generate_sigmas() actually changes the the sigmas vector inside the object.
    // Random_sigma_generator_test checks that each sigma is within the allowed range
    sigmas.set_sigmas(sigma_tot, n_sigmas, sigma_min, ebn0, esn0);
    sigmas.generate_sigmas();
    EXPECT_EQ(sigmas.get_sigmas().size(), n_sigmas) << "sigmas.generate_sigmas() should not change the size of the sigmas vector.";
    for (unsigned i = 0; i < n_sigmas; i++)
    {
        EXPECT_GE(sigmas.get_sigma(i), sigma_min) << "sigma in sigmas.get_sigmas() should be greater or equal to sigma_min.";
    }
}

TEST_F(SigmaAsymmetricTest, sigmaRatios)
{
    sigmas.set_sigmas(sigma_tot, n_sigmas, sigma_min, ebn0, esn0);
    sigmas.generate_sigmas();
    std::vector<float> sigmas_vector = sigmas.get_sigmas();
    for (unsigned i = 0; i < n_sigmas - 1; i++)
    {
        float sigma_i = sigmas_vector[i];
        float sigma_ii = sigmas_vector[i + 1];
        float expected_ratio = sigma_i * sigma_i / (sigma_i * sigma_i + sigma_ii * sigma_ii);
        EXPECT_EQ(sigmas.get_ratio(i), expected_ratio)
            << "The expected ratio of sigmas.get_ratio(i) should be the fraction between sigma of voltage level i and (i+1).";
    }
}
