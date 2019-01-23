#include <gtest/gtest.h>
#include <cmath>
#include <ql/methods/montecarlo/sample.hpp>
#include <ql/quantlib.hpp>
#include "random_number_generator_func.h"
#include "cumulative_dist_func.h"
#include "prob_density_func.h"
#include "constant.h"
#include "test_util.h"

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

/*
 * This test suite checks that the values
 * \mathbb{E}[(y-\beta_{t,y}(s))^{n}]
 * coincide between the OneSIdedBrownianBridgeFromOrigin and test_util's SignedMoment function.
 */
TEST(GenFuncOfOneSidedBrowninanBridgeTest, MonteCarloCheckForMoments){
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);
    QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> norm_rand_gen(unif_gen);

    double current_time=0.2;
    double goal_time=1.0;
    double goal_value=1.5;

    double running_sum_for_1st_moment=0;
    double running_sum_for_2nd_moment=0;
    double running_sum_for_3rd_moment=0;
    unsigned long num_of_paths=1000000;
    for(int i=0;i<num_of_paths;++i){
        double factor=0;
        double argument=0;
        random_number_generator_func::OneSidedBrownianBridgeFromOrigin(&factor,&argument,
                    goal_value,goal_time,current_time,norm_rand_gen.next().value);
        running_sum_for_1st_moment+=factor*(goal_value-argument);
        running_sum_for_2nd_moment+=factor*(goal_value-argument)*(goal_value-argument);
        running_sum_for_3rd_moment+=factor*(goal_value-argument)*(goal_value-argument)*(goal_value-argument);
    }

    EXPECT_NEAR(running_sum_for_1st_moment/num_of_paths,test_util::OneSidedBrownianBridgeSignedMoment(goal_value,goal_time,current_time,1),0.099999999);
    EXPECT_NEAR(running_sum_for_2nd_moment/num_of_paths,test_util::OneSidedBrownianBridgeSignedMoment(goal_value,goal_time,current_time,2),0.099999999);
}
