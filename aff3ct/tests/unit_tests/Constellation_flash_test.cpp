#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <string>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct;


TEST(ConstellationFlashTest, ReadIntSLCFile) {
	
	tools::Constellation_flash<int> c("test_data/test_constellation_SLC.txt");

	//TODO: Fix standard methods for completeness, not urgent.
	//EXPECT_STREQ(c.get_name(), "2Flash<C>"); 
	//EXPECT_FALSE(c.is_complex());
	
	EXPECT_EQ(c.get_n_bits_per_symbol(), 1);
	EXPECT_EQ(c.get_n_symbols(), 2);
	EXPECT_EQ(c[0], 1);
	EXPECT_EQ(c.get_symbol(1), -1);
}

TEST(ConstellationFlashTest, ReadFloatTLCFile) {
	
	tools::Constellation_flash<float> c("test_data/test_constellation_TLC.txt");

	//TODO: Fix standard methods for completeness, not urgent.
	//EXPECT_STREQ(c.get_name(), "2Flash<C>"); 
	//EXPECT_FALSE(c.is_complex());
	
	EXPECT_EQ(c.get_n_bits_per_symbol(), 3);
	EXPECT_EQ(c.get_n_symbols(), 8);
	EXPECT_EQ(c[7], 10.f);
	EXPECT_EQ(c.get_symbol(2), 6.f);
}