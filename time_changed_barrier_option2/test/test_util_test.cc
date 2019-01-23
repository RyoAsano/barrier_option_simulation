#include<gtest/gtest.h>
#include "test_util.h"
#include <ql/methods/montecarlo/sample.hpp>
#include <ql/quantlib.hpp>

TEST(TestUtilUnusualMomentTest, EqualToZero){
    EXPECT_DOUBLE_EQ(test_util::GaussianUnusualMoment(1000,0,10000,1), 0);
    EXPECT_DOUBLE_EQ(test_util::GaussianUnusualMoment(1000,10000,0,1), 0);
    EXPECT_DOUBLE_EQ(test_util::GaussianUnusualMoment(0,0,100000,2), 0);
    EXPECT_DOUBLE_EQ(test_util::GaussianUnusualMoment(0,1000000,0,2), 0);
    EXPECT_DOUBLE_EQ(test_util::GaussianUnusualMoment(1000000,0,10000000,3), 0);
    EXPECT_DOUBLE_EQ(test_util::GaussianUnusualMoment(1000000,1000000,0,3), 0);
}

TEST(TestUtilUnusualMomentTest, ExactValue){
    EXPECT_DOUBLE_EQ(test_util::GaussianUnusualMoment(1000,0.5,4.6,1), -2.3);
    EXPECT_DOUBLE_EQ(test_util::GaussianUnusualMoment(10,0.5,0.2,2), 100.01);
}

TEST(TestUtilUnusualMomentTest, MonteCarloCheck){
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);
    QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> norm_rand_gen(unif_gen);

    double rv_coeff=1.234;
    double fixed_coeff=2.33;
    double fixed_val=0.334; 
    unsigned long num_of_paths=1000000;

    double running_sum_for_1st_moment=0;
    double running_sum_for_2nd_moment=0;
    double running_sum_for_3rd_moment=0;
    double result=0;
    for(int i=0;i<num_of_paths;++i){
       double sample=norm_rand_gen.next().value*rv_coeff-fixed_coeff*fixed_val;
        running_sum_for_1st_moment+=sample;
        running_sum_for_2nd_moment+=sample*sample;
        running_sum_for_3rd_moment+=sample*sample*sample;
    }
 
    EXPECT_NEAR(running_sum_for_1st_moment/num_of_paths,test_util::GaussianUnusualMoment(rv_coeff,fixed_coeff,fixed_val,1),0.099999999);
    EXPECT_NEAR(running_sum_for_2nd_moment/num_of_paths,test_util::GaussianUnusualMoment(rv_coeff,fixed_coeff,fixed_val,2),0.099999999);
    EXPECT_NEAR(running_sum_for_3rd_moment/num_of_paths,test_util::GaussianUnusualMoment(rv_coeff,fixed_coeff,fixed_val,3),0.099999999);
}

TEST(TestUtilTruncatedMomentTest, MonteCarloCheck){
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);
    QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> norm_rand_gen(unif_gen);

    double rv_coeff=1.234;
    double fixed_coeff=2.33;
    double fixed_val=0.334; 
    unsigned long num_of_paths=1000000;

    double running_sum_for_1st_moment=0;
    double running_sum_for_2nd_moment=0;
    double running_sum_for_3rd_moment=0;
    for(int i=0;i<num_of_paths;++i){
       double sample=norm_rand_gen.next().value*rv_coeff-fixed_coeff*fixed_val;
       int integration_area_indicator=(sample>0)?1:0;
        running_sum_for_1st_moment+=sample*integration_area_indicator;
        running_sum_for_2nd_moment+=sample*sample*integration_area_indicator;
        running_sum_for_3rd_moment+=sample*sample*sample*integration_area_indicator;
    }

    EXPECT_NEAR(running_sum_for_1st_moment/num_of_paths,test_util::GaussianTruncatedMoment(rv_coeff,fixed_coeff,fixed_val,1),0.099999999);
    EXPECT_NEAR(running_sum_for_2nd_moment/num_of_paths,test_util::GaussianTruncatedMoment(rv_coeff,fixed_coeff,fixed_val,2),0.099999999);
    EXPECT_NEAR(running_sum_for_3rd_moment/num_of_paths,test_util::GaussianTruncatedMoment(rv_coeff,fixed_coeff,fixed_val,3),0.099999999);
}
   
TEST(TestUtilSignedMomentTest, MonteCarloCheck){
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);
    QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> norm_rand_gen(unif_gen);

    double rv_coeff=1.234;
    double fixed_coeff=2.33;
    double fixed_val=0.334; 
    unsigned long num_of_paths=1000000;

    double running_sum_for_1st_moment=0;
    double running_sum_for_2nd_moment=0;
    double running_sum_for_3rd_moment=0;
    for(int i=0;i<num_of_paths;++i){
       double sample=norm_rand_gen.next().value*rv_coeff-fixed_coeff*fixed_val;
       int sgn=(sample>0)?1:-1;
        running_sum_for_1st_moment+=sample*sgn;
        running_sum_for_2nd_moment+=sample*sample*sgn;
        running_sum_for_3rd_moment+=sample*sample*sample*sgn;
    }

    EXPECT_NEAR(running_sum_for_1st_moment/num_of_paths,test_util::GaussianSignedMoment(rv_coeff,fixed_coeff,fixed_val,1),0.099999999);
    EXPECT_NEAR(running_sum_for_2nd_moment/num_of_paths,test_util::GaussianSignedMoment(rv_coeff,fixed_coeff,fixed_val,2),0.099999999);
    EXPECT_NEAR(running_sum_for_3rd_moment/num_of_paths,test_util::GaussianSignedMoment(rv_coeff,fixed_coeff,fixed_val,3),0.099999999);   
}

TEST(TestUtilOneSidedBrownianBridgeSignedMoment, MonteCarloCheck){
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);
    QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> norm_rand_gen(unif_gen);

    double goal_value=1.24;
    double goal_time=2.2;
    double current_time=1.3; 
    unsigned long num_of_paths=1000000;

    double rv_coeff=sqrt((goal_time-current_time)*current_time/goal_time);
    double fixed_coeff=(goal_time-current_time)/goal_time;
    double fixed_val=goal_value;
    double factor_before_expectation=-1.0/goal_value*goal_time/(goal_time-current_time);

    double running_sum_for_1st_moment=0;
    double running_sum_for_2nd_moment=0;
    double running_sum_for_3rd_moment=0;
    for(int i=0;i<num_of_paths;++i){
       double sample=norm_rand_gen.next().value*rv_coeff-fixed_coeff*fixed_val;
       int sgn=(sample>0)?1:-1;
        running_sum_for_1st_moment+=factor_before_expectation*sample*sample*sgn;
        running_sum_for_2nd_moment+=factor_before_expectation*sample*sample*sample;
    }

    EXPECT_NEAR(running_sum_for_1st_moment/num_of_paths,test_util::OneSidedBrownianBridgeSignedMoment(goal_value,goal_time,current_time,1),0.099999999);
    EXPECT_NEAR(running_sum_for_2nd_moment/num_of_paths,test_util::OneSidedBrownianBridgeSignedMoment(goal_value,goal_time,current_time,2),0.099999999);
}

