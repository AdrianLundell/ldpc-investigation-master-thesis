#include <gtest/gtest.h>

#include <memory>
#include <vector>
#include <string>

#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

using namespace aff3ct;


TEST(ThresholderSoftTest, IntReadTestThresholds) {
	
	tools::Thresholder_soft<int> t("test_data/test_thresholder.txt");
	std::vector<int> r = {1,1,1};

	//ASSERT_EQ(t.interpret_readout(r), 1);
	//t.update_thresholds(1.f);
	//ASSERT_EQ(t.interpret_readout(r), 1);
}