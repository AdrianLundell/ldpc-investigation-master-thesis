#include <gtest/gtest.h>
#include <aff3ct.hpp>
#include <aff3ct_extension.hpp>

// Demonstrate some basic assertions.
TEST(DummyTest, BasicAssertions)
{
    // Expect two strings not to be equal.
    EXPECT_STRNE("hello", "world");
    // Expect equality.
    EXPECT_EQ(7 * 6, 42);
}
