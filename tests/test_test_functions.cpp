#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "test_functions.h"


const double eps = 1e-14;

TEST(SphereFunction_tests, exact_sol)
{
    SphereFunction func(20);
    EXPECT_NEAR(func.ymin(), func(func.solution()), eps);
}

TEST(RosenbrockFunction_tests, exact_sol)
{
    RosenbrockFunction func(20);
    EXPECT_NEAR(func.ymin(), func(func.solution()), eps);
}

TEST(RastriginFunction_tests, exact_sol)
{
    RastriginFunction func(20);
    EXPECT_NEAR(func.ymin(), func(func.solution()), eps);
}

TEST(StyblinskiTangFunction_tests, exact_sol)
{
    StyblinskiTangFunction func(20);
    EXPECT_NEAR(func.ymin(), func(func.solution()), 1e-2);
}
