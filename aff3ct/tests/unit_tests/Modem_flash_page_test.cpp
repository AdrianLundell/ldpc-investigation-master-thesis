#include <gtest/gtest.h>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

namespace aff3ct 
{

TEST(ModemFlashPageTest, slc)
{
    int N = 3;
    tools::Flash_cell cell = tools::Flash_cell(Flash_cell::SLC);
    tools::Flash_reader reader = tools::Flash_reader(tools::Flash_reader::lower, tools::Flash_reader::hard);
    tools::Noise<R>& noise;

    module::Modem_flash_page m = module::Modem_flash_page(N, cell, reader, noise);

}




}