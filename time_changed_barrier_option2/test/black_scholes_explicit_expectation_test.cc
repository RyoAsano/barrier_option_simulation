#include <gtest/gtest.h>
#include "black_scholes_explicit_expectation.h"
#include <cmath>

TEST(BlackScholesExplicitExpectationTest, ExactValue){
    double drift=0.1;
    double vol=0.2;
    double maturity=0.5;
    double strike=40;
    double init_val=42;

    EXPECT_NEAR(exp(-0.05)*black_scholes_explicit_expectation::VanillaCall(
                drift,vol,maturity,strike,init_val), 4.76, 0.00999999999);

    EXPECT_NEAR(exp(-0.05)*black_scholes_explicit_expectation::VanillaPut(
                drift,vol,maturity,strike,init_val), 0.81, 0.00999999999);
}

TEST(BlackScholesUpAndInExpectationTest, monotone_decreasing_wrt_barrier_level){
    bool monotone_decreasing_wrt_barrier_level=true;
    bool monotone_increasing_wrt_initial_value=true;
    bool monotone_decreasing_wrt_strike=true;
    bool monotone_increasing_wrt_drfit=true;

    double drift_init=-0.2;
    double volatility=0.2;
    double maturity=1.0;
    double strike_init=100;
    double initial_value=110;
    double barrier_level_init=initial_value;

    double drift_max=0.2;
    double barrier_level_max=120;
    double strike_max=barrier_level_max;
    double num_of_subdivisions=10;

    double drift_running=drift_init;
    for(int i=0;i<num_of_subdivisions;++i){
        double strike_running=strike_init;
        double barrier_level_running=barrier_level_init;
        for(int j=0;j<num_of_subdivisions;++j){
            barrier_level_running=barrier_level_init;    
            for(int k=0;k<num_of_subdivisions;++k){
                double before=black_scholes_explicit_expectation::UpAndInCall(drift_running,
                        volatility,maturity,strike_running,initial_value,barrier_level_running);
                double barrier_level_before=barrier_level_running;
                barrier_level_running+=(barrier_level_max-barrier_level_init)/num_of_subdivisions;
                double after=black_scholes_explicit_expectation::UpAndInCall(drift_running,
                        volatility,maturity,strike_running,initial_value,barrier_level_running);
                if(strike_running<barrier_level_before){
                    EXPECT_GE(before,after)
                        << "barrier_level:"<<barrier_level_before<<"->"<<barrier_level_running
                        <<" drift:"<<drift_running<<" strike:"<<strike_running; 
                }else{
                    EXPECT_NEAR(before,after,0.0099999999)
                        << "barrier_level:"<<barrier_level_before<<"->"<<barrier_level_running
                        <<" drift:"<<drift_running<<" strike:"<<strike_running;
                }
           }
            strike_running+=(strike_max-strike_init)/num_of_subdivisions;
        }
        drift_running+=(drift_max-drift_init)/num_of_subdivisions; 
    }
}

TEST(BlackScholesUpAndInExpectationTest, monotone_increasing_wrt_drfit){
    double drift_init=-0.2;
    double volatility=0.2;
    double maturity=1.0;
    double strike_init=100;
    double initial_value=110;
    double barrier_level_init=initial_value;

    double drift_max=0.2;
    double barrier_level_max=120;
    double strike_max=barrier_level_max;
    double num_of_subdivisions=10;

    double barrier_level_running=barrier_level_init;
    for(int i=0;i<num_of_subdivisions;++i){
        double strike_running=strike_init;
        for(int j=0;j<num_of_subdivisions;++j){
            double drift_running=drift_init;
            for(int k=0;k<num_of_subdivisions;++k){
                double before=black_scholes_explicit_expectation::UpAndInCall(drift_running,
                        volatility,maturity,strike_running,initial_value,barrier_level_running);
                double drift_running_before=drift_running;
                drift_running+=(drift_max-drift_init)/num_of_subdivisions;
                double after=black_scholes_explicit_expectation::UpAndInCall(drift_running,
                        volatility,maturity,strike_running,initial_value,barrier_level_running);
                EXPECT_LE(before,after)
                    <<"drift:"<<drift_running_before<<"->"<<drift_running
                    <<" strike:"<<strike_running<<" barrier_level:"<<barrier_level_running; 
            }
            strike_running+=(strike_max-strike_init)/num_of_subdivisions;
        }
        barrier_level_running+=(barrier_level_max-barrier_level_init)/num_of_subdivisions;
    }
}

TEST(BlackScholesUpAndInExpectationTest, monotone_decreasing_wrt_strike){
    double drift_init=-0.2;
    double volatility=0.2;
    double maturity=1.0;
    double strike_init=100;
    double initial_value=110;
    double barrier_level_init=initial_value;

    double drift_max=0.2;
    double barrier_level_max=120;
    double strike_max=barrier_level_max;
    double num_of_subdivisions=10;

    double drift_running=drift_init;
    for(int i=0;i<num_of_subdivisions;++i){
        double barrier_level_running=barrier_level_init;
        for(int j=0;j<num_of_subdivisions;++j){
            double strike_running=strike_init;
            for(int k=0;k<num_of_subdivisions;++k){
                double before=black_scholes_explicit_expectation::UpAndInCall(drift_running,
                        volatility,maturity,strike_running,initial_value,barrier_level_running);
                double strike_running_before=strike_running;
                strike_running+=(strike_max-strike_init)/num_of_subdivisions;
                double after=black_scholes_explicit_expectation::UpAndInCall(drift_running,
                        volatility,maturity,strike_running,initial_value,barrier_level_running);
                EXPECT_GE(before,after)
                    <<"strike:"<<strike_running_before<<"->"<<strike_running
                    <<" barrier_level:"<<barrier_level_running<<" drift:"<<drift_running;
           }
            barrier_level_running+=(barrier_level_max-barrier_level_init)/num_of_subdivisions;
        }
        drift_running+=(drift_max-drift_init) /num_of_subdivisions;
    }
}

TEST(BlackScholesUpAndInExpectationTest,should_equal_vanilla_call){
    double drift=0.02;
    double volatility=0.2;
    double maturity=1.0;
    double strike=100;
    double initial_value=100;
    double barrier_level=initial_value;

    double vanilla=black_scholes_explicit_expectation::VanillaCall(drift,volatility,maturity,strike,initial_value);
    double up_and_in=black_scholes_explicit_expectation::UpAndInCall(drift,volatility,maturity,strike,initial_value,barrier_level);

    EXPECT_DOUBLE_EQ(vanilla,up_and_in);
}

TEST(BlackScholesCallDeltaTest, Benchmark){
    double drift=0;
    double volatility=0.2;
    double maturity=0.3;
    double strike=100;
    double initial_value=100;

    for(int i=0;i<10;++i){
        std::cerr << "[          ] volatility:"<< volatility 
                  <<" delta:" << black_scholes_explicit_expectation::VanillaCallDelta(drift,volatility,maturity,strike, initial_value)
                  << std::endl; 
        volatility+=0.05;
    }
}
