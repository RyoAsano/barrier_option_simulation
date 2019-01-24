#include <gtest/gtest.h>
#include "prob_density_func.h"
#include <ql/methods/montecarlo/sample.hpp>
#include <ql/quantlib.hpp>

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

TEST(ProbDensityOfStdGaussian, HalfRealLineMeasureShouldEqualOneHalf){
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);

    double running_sum_for_integration=0;
    unsigned long num_of_paths=10000000;
    for(int i=0;i<num_of_paths;++i){
        double unif=unif_gen.next().value;
        double changed_variable=0-(1.0-unif)/unif;
        double integrator_multiplier=1.0/(unif*unif);
        running_sum_for_integration+=prob_density_func::StdGaussian(changed_variable)*integrator_multiplier;
    }
    EXPECT_NEAR(running_sum_for_integration/num_of_paths,0.5,0.0999999);
}

TEST(ProbDensityOfGaussian, ExactValue){ 
    EXPECT_NEAR(prob_density_func::Gaussian(0.234, 0.586*0.586, 1.23),0.16058644,0.00000000999999999);
}

TEST(ProbDensityOfGaussian, HalfRealLineMeasureShouldEqualOneHalf){
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);

    double mean=1.2;
    double variance=2.4;
    double running_sum_for_integration=0;
    unsigned long num_of_paths=10000000;
    for(int i=0;i<num_of_paths;++i){
        double unif=unif_gen.next().value;
        double changed_variable=mean-(1.0-unif)/unif;
        double integrator_multiplier=1.0/(unif*unif);
        running_sum_for_integration+=prob_density_func::Gaussian(mean,variance,changed_variable)*integrator_multiplier;
    }
    EXPECT_NEAR(running_sum_for_integration/num_of_paths,0.5,0.0999999);
}

TEST(ProbDensistyOfOneSidedBB, TotalMeasureShouldEqualOne){

    double goal_value=1.5;
    double goal_time=0.7;
    double current_time=0.3845;

    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);

    double running_sum_for_integration=0;
    unsigned long num_of_paths=10000000;
    for(int i=0;i<num_of_paths;++i){
        double unif=unif_gen.next().value;
        double changed_variable=goal_value-(1.0-unif)/unif;
        double integrator_multiplier=1.0/(unif*unif);
        running_sum_for_integration+=prob_density_func::OneSidedBrownianBridge(goal_value,goal_time,
                current_time,changed_variable)*integrator_multiplier;
    }
    EXPECT_NEAR(running_sum_for_integration/num_of_paths,1.0,0.0999999);
}

