#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <string>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct;


TEST(ThresholderSoftTest, IntReadTestThresholds) {
	
	tools::Thresholder_soft<int> t("test_data/test_thresholder.txt");
	t.update_thresholds();

	//Int -> rounded interpolation
	ASSERT_EQ(t.get_threshold(0), 1);
	ASSERT_EQ(t.get_threshold(1), 3);
	ASSERT_EQ(t[2], 5);

	std::vector<int> r = {1,1,1};
	ASSERT_EQ(t.interpret_readout(r), 15);
	
	r = {0,1,1};
	ASSERT_EQ(t.interpret_readout(r), 35);
	
	r = {0,0,1};
	ASSERT_EQ(t.interpret_readout(r), 55);
	
	r = {0,0,0};
	ASSERT_EQ(t.interpret_readout(r), 75);
}


TEST(ThresholderSoftTest, FloatReadTestThresholds) {
	
	tools::Thresholder_soft<float> t("test_data/test_thresholder.txt");
	t.update_thresholds();

	//Int -> rounded interpolation
	ASSERT_FLOAT_EQ(t.get_threshold(0), 1.5f);
	ASSERT_FLOAT_EQ(t.get_threshold(1), 3.5f);
	ASSERT_FLOAT_EQ(t[2], 5.5f);

	std::vector<int> r = {1,1,1};
	ASSERT_FLOAT_EQ(t.interpret_readout(r), 15.f);
	
	r = {0,1,1};
	ASSERT_FLOAT_EQ(t.interpret_readout(r), 35.f);
	
	r = {0,0,1};
	ASSERT_FLOAT_EQ(t.interpret_readout(r), 55.f);
	
	r = {0,0,0};
	ASSERT_FLOAT_EQ(t.interpret_readout(r), 75.f);
}