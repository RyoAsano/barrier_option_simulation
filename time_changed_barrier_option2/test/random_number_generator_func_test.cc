#include <gtest/gtest.h>
#include <cmath>
#include "random_number_generator_func.h"
#include "constant.h"

TEST(IverseCDFofExpDistTest, ExactValue){
    EXPECT_NEAR(random_number_generator_func::Exponential(1.0, 0.2), 1.609437912, 0.00000000099999999);
    EXPECT_NEAR(random_number_generator_func::Exponential(3.0, 0.46), 0.2588429298, 0.00000000099999999);
}

/*
 * note that the downward cdf of exp dist. (i.e. P[E>x]) is exp(-lambda*x)
 */
TEST(InverseCDFofExpDistTest, CheckInverseRelation){     
    EXPECT_DOUBLE_EQ(random_number_generator_func::Exponential(1.0, exp(-1.0*3.0)), 3.0);
    EXPECT_DOUBLE_EQ(random_number_generator_func::Exponential(1.0, exp(-1.0*6.0)), 6.0);
    EXPECT_DOUBLE_EQ(random_number_generator_func::Exponential(3.0, exp(-3.0*2.3)), 2.3);
}

TEST(InverseCDFofArcsineDistTest, CheckInverseRelation){ 
    EXPECT_NEAR(random_number_generator_func::Arcsine(1.0-2.0/PI*asin(sqrt(0.4))), 0.4, 0.000000001);
    EXPECT_NEAR(random_number_generator_func::Arcsine(1.0-2.0/PI*asin(sqrt(0.23))), 0.23, 0.0000000001);
}

TEST(GenFuncOfGammaShapeAHalfScaleOneTest, CheckInverseRelation){ 
    double unif1=1.0-2.0/PI*asin(sqrt(0.4));
    double unif2=exp(-1.0*3.4);
    EXPECT_DOUBLE_EQ(random_number_generator_func::GammaShapeAHalfScaleOne(unif1, unif2), 0.4*3.4);
}

TEST(GenFuncOfBMsFirstHittingTimeTest, CheckInverseRelation){   
    double a = 0.4;
    double b = 3.4;
    double barrier_level=3.0;
    double unif1=1.0-2.0/PI*asin(sqrt(a));
    double unif2=exp(-1.0*b);
    double expected_value = barrier_level*barrier_level/(2.0*a*b);
    EXPECT_DOUBLE_EQ(random_number_generator_func::BrownianMotionFirstHittingTime(barrier_level,unif1,unif2),expected_value);
}
