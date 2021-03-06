#include <gtest/gtest.h>
#include "cumulative_dist_func.h"

TEST(CDFofStdGaussianDist, ProbHalf){ 
    EXPECT_DOUBLE_EQ(0.5, cumulative_dist_func::StdGaussian(0));
    EXPECT_LT(cumulative_dist_func::StdGaussian(-0.05), 0.5);
    EXPECT_LT(0.5, cumulative_dist_func::StdGaussian(0.05));
}

TEST(CDFofStdGaussianDist, Symmetry){
    EXPECT_DOUBLE_EQ(cumulative_dist_func::StdGaussian(0.5), 1.0-cumulative_dist_func::StdGaussian(-0.5));
    EXPECT_DOUBLE_EQ(cumulative_dist_func::StdGaussian(1.0), 1.0-cumulative_dist_func::StdGaussian(-1.0));
    EXPECT_DOUBLE_EQ(cumulative_dist_func::StdGaussian(-0.2), 1.0-cumulative_dist_func::StdGaussian(0.2));
}

TEST(CDFofStdGaussianDist, ExactValue){
    EXPECT_NEAR(cumulative_dist_func::StdGaussian(-2.82), 0.0024, 0.0000999999999); 
    EXPECT_NEAR(cumulative_dist_func::StdGaussian(-1.22), 0.1112, 0.0000999999999);
    EXPECT_NEAR(cumulative_dist_func::StdGaussian(1.58), 0.9429, 0.0000999999999);
    EXPECT_NEAR(cumulative_dist_func::StdGaussian(0.45), 0.6736, 0.0000999999999); 
}
