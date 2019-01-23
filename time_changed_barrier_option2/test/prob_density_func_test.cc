#include <gtest/gtest.h>
#include "prob_density_func.h"

TEST(ProbDensityOfBrownainMotionFirstHittingTimeTest, TimeZero){
    EXPECT_DOUBLE_EQ(0, prob_density_func::BrownianMotionFirstHittingTime(1.0, 0));
    EXPECT_DOUBLE_EQ(0, prob_density_func::BrownianMotionFirstHittingTime(0,0));
}

TEST(ProbDensityOfBrownainMotionFirstHittingTimeTest, ExactValue){
    EXPECT_NEAR(prob_density_func::BrownianMotionFirstHittingTime(3.0,1.5), 0.03243478222, 0.00000000000999999999);
    EXPECT_NEAR(prob_density_func::BrownianMotionFirstHittingTime(2.2,9.5), 0.02323360138, 0.00000000000999999999);
}

TEST(ProbDensityOfStdGaussian, ExactValue){ 
    EXPECT_NEAR(prob_density_func::StdGaussian(0.323),0.37866513,0.00000000999999999);
}

