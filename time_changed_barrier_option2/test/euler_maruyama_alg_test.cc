#include <gtest/gtest.h>
#include "euler_maruyama_alg.h"
#include "test_util.h" 

TEST(TimeChangedBlackScholesUpperBarrierOneStepTest, DeterministicLinear){
    double drift=0.02;
    double volatility=0.2;
    double asset_price_process_init=1.0;
    double asset_price_process=asset_price_process_init;
    double barrier_hitting_time=0;
    double goal_value=2.0;
    double goal_time=2.0;


    double num_of_subdivisions=100000;
    double time_increment=goal_time/num_of_subdivisions;
    double deterministic_linear_integrator_increment=(goal_value-asset_price_process_init)*time_increment/goal_time;
    double current_time=0;
    for(int i=0;i<num_of_subdivisions;++i){
        euler_maruyama::TimeChangedBlackScholesUpperBarrierOneStep(volatility,deterministic_linear_integrator_increment,
                time_increment,&asset_price_process,&barrier_hitting_time);
        double expected_asset_price=0;
        double expected_barrier_hitting_time=0;
        test_util::BlackScholesTimeChangedHittingTimeLinearCase(volatility,asset_price_process_init,goal_value,
                goal_time,current_time,&expected_asset_price,&expected_barrier_hitting_time);
        EXPECT_NEAR(asset_price_process,expected_asset_price,0.0009999);
        EXPECT_NEAR(barrier_hitting_time,expected_barrier_hitting_time,0.000999);
        current_time+=time_increment;
    }
}
