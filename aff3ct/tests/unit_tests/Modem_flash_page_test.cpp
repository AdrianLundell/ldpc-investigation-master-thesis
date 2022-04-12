#include <gtest/gtest.h>
#include <aff3ct.hpp> 
#include <aff3ct_extension.hpp>

namespace aff3ct
{

TEST(ModemFlashPageTest, slc)
{
    tools::Flash_cell cell = Flash_cell(tools::Flash_cell::SLC);
    //Set symmetric noise
    tools::Flash_reader<float, float> reader(tools::Flash_reader<float, float>::lower, tools::Flash_reader<float, float>::hard, "test_data/reader_static_hard.txt");
    tools::Sigma_asymmetric<float> s1;
    std::vector<float> sigmas1 = {1, 1};
    s1.set_sigmas(sigmas1);
    std::vector<float> levels1 = {-1, 1};
    aff3ct::module::Channel_AWGN_asymmetric<float> c1(1, levels1, s1, 0, 1);
    c1.set_noise(s1);
    reader.update(c1);

    module::Modem_flash_page<float, float, float> f(2, cell, reader, s1, 1);    

    std::vector<float> codeword = {0,1};
    std::vector<float> symbols = {0, 0};
    f.modulate(codeword, symbols);

    ASSERT_FLOAT_EQ(symbols[0], 0);
    ASSERT_FLOAT_EQ(symbols[1], 1);

    std::vector<float> noisy_symbols = {1.1, -1.1};
    std::vector<float> output = {0, 1};
    f.demodulate(noisy_symbols, output);

    ASSERT_FLOAT_EQ(symbols[0], 0);
    ASSERT_FLOAT_EQ(symbols[1], 1);
}
}