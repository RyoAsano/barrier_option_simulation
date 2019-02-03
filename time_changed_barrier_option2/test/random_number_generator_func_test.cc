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

TEST(GenFuncOfBMsFirstHittingTimeTest, MonteCarloCheck){
    double maturity=1.0;
    double barrier_level=1.0;
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed); 
    QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> norm_rand_gen(unif_gen);

    double running_sum_for_euler_maruyama=0;
    double running_sum_for_first_hitting_time=0;

    unsigned long num_of_paths=10000;
    unsigned int num_of_subdivisions=10;
    double time_increment=maturity/num_of_subdivisions;
    for(int i=0;i<num_of_paths;++i){
        double running_brownian_motion=0;
        double indicator_function=0;
        for(int k=0;k<num_of_subdivisions;++k){
            if(running_brownian_motion>barrier_level){
              indicator_function=1.0;  
              break;
            }
            running_brownian_motion+=norm_rand_gen.next().value*sqrt(time_increment);
        }
        running_sum_for_euler_maruyama+=indicator_function; 
    }

    for(int i=0;i<num_of_paths;++i){
        double sample_first_hitting_time=random_number_generator_func::BrownianMotionFirstHittingTime(barrier_level,
                unif_gen.next().value,unif_gen.next().value);
        if(sample_first_hitting_time<=maturity){
            running_sum_for_first_hitting_time+=1.0;
        }
    }
    EXPECT_NEAR(running_sum_for_first_hitting_time/num_of_paths,running_sum_for_euler_maruyama/num_of_paths,0.099999);
}

/*
 * This test suite checks that the values
 * \mathbb{E}[(y-\beta_{t,y}(s))^{n}]
 * coincide between the OneSIdedBrownianBridgeFromOrigin and test_util's Moment function.
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

    EXPECT_NEAR(running_sum_for_1st_moment/num_of_paths,test_util::OneSidedBrownianBridgeMoment(goal_value,goal_time,current_time,1),0.099999999);
    EXPECT_NEAR(running_sum_for_2nd_moment/num_of_paths,test_util::OneSidedBrownianBridgeMoment(goal_value,goal_time,current_time,2),0.099999999);
}

TEST(GenFuncOfOneSidedBrowninanBridgeTest,IncrementShouldEqualFromOrigin){
    double goal_value=1.872;
    double goal_time=7.237;
    double std_norm_rand_num=0.3333;

    double current_time_for_from_origin=0.118;
    double factor_for_from_origin=0;
    double arg_for_from_origin=0; 
    random_number_generator_func::OneSidedBrownianBridgeFromOrigin(&factor_for_from_origin,&arg_for_from_origin,
            goal_value,goal_time,current_time_for_from_origin,std_norm_rand_num);

    double current_time_for_increment=0;
    double time_increment=current_time_for_from_origin-current_time_for_increment;
    double factor_for_increment=1; //Note that in case of use of increment function, you must initialize it with 1.
    double arg_for_increment=0;
    arg_for_increment=random_number_generator_func::OneSidedBrownianBridgeIncrement(goal_value,goal_time,
            current_time_for_increment,0,time_increment,std_norm_rand_num, &factor_for_increment);
    EXPECT_DOUBLE_EQ(factor_for_from_origin,factor_for_increment);
    EXPECT_DOUBLE_EQ(arg_for_from_origin,arg_for_increment);
}

TEST(GenFuncOfOneSidedBrowninanBridgeTest, UsingDensityFunction){
    double goal_value=1.434;
    double goal_time=0.921;
    double current_time=0.223;
    auto test_function =[](double argument){return (argument>0)?argument:0;};

    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);
    QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> norm_rand_gen(unif_gen);

    double running_sum_for_generator=0;
    double running_sum_for_density=0;
    unsigned long num_of_paths=1000000;
    for(int i=0;i<num_of_paths;++i){
        double factor=1.0;
        double argument=0;
        random_number_generator_func::OneSidedBrownianBridgeFromOrigin(&factor,&argument,
                goal_value,goal_time,current_time,norm_rand_gen.next().value);
        running_sum_for_generator+=factor*test_function(argument);
        
        double unif=unif_gen.next().value;
        double changed_variable=goal_value-(1.0-unif)/unif;
        double integrator_multiplier=1.0/(unif*unif);
        running_sum_for_density+=test_function(changed_variable)*
            prob_density_func::OneSidedBrownianBridge(goal_value,goal_time,
                current_time,changed_variable)*integrator_multiplier;
    }
    EXPECT_NEAR(running_sum_for_density/num_of_paths,running_sum_for_generator/num_of_paths,0.099999);
}
