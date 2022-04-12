#include <gtest/gtest.h>
#include <aff3ct_extension.hpp>
#include <aff3ct.hpp>

#include <iostream>
#include <vector>

using namespace aff3ct;
using std::vector;

class ChannelAwgnAsymmetricTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        N = 2;
        ebn0 = 0.01f;
        R = 1;
        esn0 = tools::ebn0_to_esn0(ebn0, R);
        sigma_tot = tools::esn0_to_sigma(esn0);
        sigma_min = 0.05;
        voltage_levels = {-5.0f, 5.0f};
        n_sigmas = voltage_levels.size();
        sigmas.set_sigmas(sigma_tot, n_sigmas, sigma_min, ebn0, esn0);
        sigmas.set_seed(3);
        sigmas.generate_sigmas();
        voltage_level_indexes = {0, 0};
        noisy_voltages = std::vector<float>(N);
    }

    int N;
    float ebn0;
    float R;
    float esn0;
    float sigma_tot;
    float sigma_min;
    unsigned n_sigmas;
    vector<float> voltage_levels;
    tools::Sigma_asymmetric<float> sigmas;
    vector<unsigned> voltage_level_indexes;
    vector<float> noisy_voltages;
};

TEST_F(ChannelAwgnAsymmetricTest, initialize)
{
    module::Channel_AWGN_asymmetric<float> channel(N, voltage_levels, sigmas);
    EXPECT_NO_THROW(channel.check_noise()) << "Channel noise should be of type SIGMA.";
}

TEST_F(ChannelAwgnAsymmetricTest, addNoise)
{
    module::Channel_AWGN_asymmetric<float> channel(N, voltage_levels, sigmas);
    channel.add_noise(voltage_level_indexes.data(), noisy_voltages.data());
    for (float noisy_voltage_level : noisy_voltages)
    {
        std::cout << noisy_voltage_level << std::endl;
        EXPECT_LE(noisy_voltage_level, 0.0f);
    }
}