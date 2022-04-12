#include <gtest/gtest.h>
#include <aff3ct_extension.hpp>
#include <aff3ct.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <random>

using namespace aff3ct;
using std::vector;

class ChannelAwgnAsymmetricTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        printData = true; // Set to true to print randomly generated sigmas and noise
        N = 4;
        ebn0 = 0.01f;
        R = 1;
        esn0 = tools::ebn0_to_esn0(ebn0, R);
        sigma_tot = tools::esn0_to_sigma(esn0);
        sigma_min = 0.05;
        voltage_levels = {-5.0f, 5.0f, 10.0f}; // First value should be well below zero for addNoise to succeed
        n_sigmas = voltage_levels.size();
        sigmas.set_sigmas(sigma_tot, n_sigmas, sigma_min, ebn0, esn0);
        sigmas.generate_sigmas();
        voltage_level_indexes = {0, 0, 0, 0}; // Should be all zeros for addNoise to succeed
        noisy_voltages = std::vector<float>(N);
    }

    void writeData(module::Channel_AWGN_asymmetric<float> &channel)
    {
        // Make multiple runs of add_noise and save the sigmas and the noise in separate text files
        std::ofstream sigma_file("../data/sigmas.txt");
        std::ofstream noise_file("../data/noise.txt");
        std::ofstream parameters_file("../data/parameters.txt");
        if (sigma_file.is_open() && noise_file && parameters_file)
        {
            parameters_file << N << "\n"
                            << sigma_tot << "\n"
                            << sigma_min << "\n"
                            << n_sigmas << "\n"
                            << voltage_levels[0] << "\n"
                            << sigmas.get_sigma((unsigned)0) << "\n";

            for (unsigned n = 0; n < 10000; n++)
            {
                channel.add_noise(voltage_level_indexes.data(), noisy_voltages.data());

                for (float &noisy_voltage : noisy_voltages)
                {
                    noise_file << noisy_voltage << "\n";
                }
            }

            for (unsigned n = 0; n < 100000; n++)
            {
                channel.generate_sigmas();
                for (unsigned i = 0; i < n_sigmas; i++)
                {
                    sigma_file << channel.get_sigma(i)
                               << "\n";
                }

                for (float &noisy_voltage : noisy_voltages)
                {
                    noise_file << noisy_voltage << "\n";
                }
            }
        }
        sigma_file.close();
        noise_file.close();
        parameters_file.close();
    }

    bool printData;
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
        EXPECT_LE(noisy_voltage_level, 0.0f);
    }
    if (printData)
    {
        this->writeData(channel);
    }
}