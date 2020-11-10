#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "test_functions.h"


const double eps = 1e-14;

TEST(SphereFunction_tests, exact_sol)
{
    SphereFunction<20> func;
    EXPECT_NEAR(func.ymin(), func(func.solution())[0], eps);
}

TEST(RosenbrockFunction_tests, exact_sol)
{
    RosenbrockFunction<20> func;
    EXPECT_NEAR(func.ymin(), func(func.solution())[0], eps);
}

TEST(RastriginFunction_tests, exact_sol)
{
    RastriginFunction<20> func;
    EXPECT_NEAR(func.ymin(), func(func.solution())[0], eps);
}

TEST(StyblinskiTangFunction_tests, exact_sol)
{
    StyblinskiTangFunction<20> func;
    EXPECT_NEAR(func.ymin(), func(func.solution())[0], 1e-2);
}
