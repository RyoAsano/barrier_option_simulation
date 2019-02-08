#include <gtest/gtest.h>
#include "euler_maruyama_alg.h"
#include "test_util.h" 
#include <ql/methods/montecarlo/sample.hpp>
#include <ql/quantlib.hpp>
#include "black_scholes_explicit_expectation.h" 
#include <unistd.h>
#include <time.h>
#include "random_number_generator_func.h"
#include "statistical_functions.h"


TEST(TimeChangedBlackScholesUpperBarrierOneStepTest, DeterministicLinear){
    double drift=0;
    double volatility=0.2;
    double asset_price_process_init=1.0;
    double asset_price_process=asset_price_process_init;
    double barrier_hitting_time=0;
    double log_density_process=0;
    double goal_value=2.0;
    double goal_time=2.0;

    double num_of_subdivisions=100000;
    double time_increment=goal_time/num_of_subdivisions;
    double deterministic_linear_integrator_increment=(goal_value-asset_price_process_init)*time_increment/goal_time;
    double current_time=0;
    for(int i=0;i<num_of_subdivisions;++i){
        euler_maruyama::TimeChangedBlackScholesUpperBarrierOneStep(drift,volatility,deterministic_linear_integrator_increment,
                time_increment,&asset_price_process,&barrier_hitting_time,&log_density_process);
        double expected_asset_price=0;
        double expected_barrier_hitting_time=0;
        test_util::BlackScholesTimeChangedHittingTimeLinearCase(volatility,asset_price_process_init,goal_value,
                goal_time,current_time,&expected_asset_price,&expected_barrier_hitting_time);
        EXPECT_NEAR(asset_price_process,expected_asset_price,0.0009999);
        EXPECT_NEAR(barrier_hitting_time,expected_barrier_hitting_time,0.000999);
        current_time+=time_increment;
    }
}

TEST(test,test){
    double array[10000]={0};
    double num=10000;
    EXPECT_EQ(sizeof(array)/sizeof(array[0]),num);
}

TEST(TimeChangedBlackScholesUpperBarrierOnePathTest,DeterministicLinear){
    double drift=0;
    double volatility=0.2;
    double maturity=1000000;
    double sde_initial_value=100;
    double barrier_level=120;

    double asset_price_process;
    double barrier_hitting_time;
    double log_density_process;
    double factor;
    bool the_process_hits_the_barrier_before_maturity;

    double expected_asset_price;
    double expected_barrier_hitting_time;
    double brownian_bridge_goal_value=barrier_level;

    double norm_rand_num_array[100]={0};
    unsigned long num_of_subdivisions=sizeof(norm_rand_num_array)/sizeof(norm_rand_num_array[0]);
    double running_goal_time=0;
    double brownian_bridge_goal_time=2.0;
    double num_of_subdivisions_for_running_goal_time=100000;
    double time_increment=brownian_bridge_goal_time/num_of_subdivisions_for_running_goal_time;
    assert((sizeof(norm_rand_num_array)/sizeof(norm_rand_num_array[0]))==num_of_subdivisions);

    for(int i=0;i<num_of_subdivisions_for_running_goal_time;++i){
        running_goal_time+=time_increment;
        euler_maruyama::TimeChangedBlackScholesUpperBarrierOnePath(drift,volatility,maturity,sde_initial_value,barrier_level,
                running_goal_time,num_of_subdivisions,norm_rand_num_array,&asset_price_process,&barrier_hitting_time,&log_density_process,
                &factor,&the_process_hits_the_barrier_before_maturity);

        test_util::BlackScholesTimeChangedHittingTimeLinearCase(volatility,sde_initial_value,brownian_bridge_goal_value,
                running_goal_time,running_goal_time,&expected_asset_price,&expected_barrier_hitting_time);

        //reason why there is a certain discrepancy is that the OnePath function uses the integration approximation for the hitting time.
        EXPECT_NEAR(barrier_hitting_time,expected_barrier_hitting_time,0.0099999);
    }
}

TEST(TimeChangedBlackScholesUpperBarrierOnePathTest, BeTheGoalValueAtTheGoalTime){
    double drift=0;
    double volatility=0.2;
    double maturity=100000;
    double sde_initial_value=100;
    double barrier_level=120;
    double brownian_bridge_goal_time=1.0;

    double asset_price_process;
    double barrier_hitting_time;
    double log_density_process;
    double factor;
    bool the_process_hits_the_barrier_before_maturity;

    unsigned long num_of_subdivisions=1000;
    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);
    QuantLib::RandomSequenceGenerator<QuantLib::MersenneTwisterUniformRng> norm_rand_num_array(num_of_subdivisions,unif_gen);

    euler_maruyama::TimeChangedBlackScholesUpperBarrierOnePath(drift,volatility,maturity,sde_initial_value,barrier_level,
            brownian_bridge_goal_time,num_of_subdivisions,&(norm_rand_num_array.nextSequence().value[0]),
            &asset_price_process,&barrier_hitting_time,&log_density_process,&factor,&the_process_hits_the_barrier_before_maturity);
    EXPECT_NEAR(asset_price_process,barrier_level,0.0000000099999);
}

TEST(TimeChangedBlackScholesUpperBarrierOnePathTest, CheckNotHittingTheBarrierInDeterministicLinearCase){
    double drift=0;
    double volatility=0.2;
    double sde_initial_value=100;
    double barrier_level=120;

    double asset_price_process;
    double barrier_hitting_time;
    double log_density_process;
    double factor;
    bool the_process_hits_the_barrier_before_maturity;

    double expected_asset_price;
    double expected_barrier_hitting_time;
    double brownian_bridge_goal_value=barrier_level;

    double norm_rand_num_array[100]={0};
    unsigned long num_of_subdivisions=sizeof(norm_rand_num_array)/sizeof(norm_rand_num_array[0]);
    double running_goal_time=0;
    double brownian_bridge_goal_time=2.0;
    double num_of_subdivisions_for_running_goal_time=100000;
    double time_increment=brownian_bridge_goal_time/num_of_subdivisions_for_running_goal_time;
    assert((sizeof(norm_rand_num_array)/sizeof(norm_rand_num_array[0]))==num_of_subdivisions);

    for(int i=0;i<num_of_subdivisions_for_running_goal_time;++i){
        running_goal_time+=time_increment;
        test_util::BlackScholesTimeChangedHittingTimeLinearCase(volatility,sde_initial_value,brownian_bridge_goal_value,
                running_goal_time,running_goal_time,&expected_asset_price,&expected_barrier_hitting_time);
        double maturity=expected_barrier_hitting_time-0.0099999999;

        euler_maruyama::TimeChangedBlackScholesUpperBarrierOnePath(drift,volatility,maturity,sde_initial_value,barrier_level,
                running_goal_time,num_of_subdivisions,norm_rand_num_array,&asset_price_process,&barrier_hitting_time,&log_density_process,
                &factor,&the_process_hits_the_barrier_before_maturity);
        EXPECT_FALSE(the_process_hits_the_barrier_before_maturity)
            <<"estimated hitting time:"<<barrier_hitting_time<<" true hitting time:"<<expected_barrier_hitting_time;
    }
}

TEST(TimeChangedBlackScholesUpperBarrierOnePathTest,MonteCarloCheckForSimpleExpectation){
    double drift=0;
    double volatility=0.2;
    double sde_initial_value=100;
    double barrier_level=101;
    double maturity=1.0;
    unsigned int  num_of_subdivisions=10;
    unsigned long num_of_paths=100000;

    double asset_price_process;
    double barrier_hitting_time;
    double factor;
    double the_process_hits_the_barrier_before_maturity;

    time_t now=time(0);
    unsigned long seed=getpid()+(unsigned long)now;

    typedef QuantLib::MersenneTwisterUniformRng MersenneTwister;
    typedef QuantLib::BoxMullerGaussianRng<MersenneTwister> BoxMuller;
    QuantLib::RandomSequenceGenerator<MersenneTwister> unif_gen(2,MersenneTwister(seed));
    QuantLib::RandomSequenceGenerator<BoxMuller> norm_rand_num_array_gen(num_of_subdivisions,BoxMuller(MersenneTwister(seed)));

    double brownian_bridge_goal_time=1.0;

    double running_sum_for_monte_carlo=0;
    for(int i=0;i<num_of_paths;++i){
        double asset_price_process;
        double barrier_hitting_time;
        double log_density_process;
        double factor;
        bool the_process_hits_the_barrier_before_maturity;
        euler_maruyama::TimeChangedBlackScholesUpperBarrierOnePath(drift,volatility,maturity,sde_initial_value,
                barrier_level,brownian_bridge_goal_time,num_of_subdivisions,&(norm_rand_num_array_gen.nextSequence().value[0]),
                &asset_price_process,&barrier_hitting_time,&log_density_process,&factor,&the_process_hits_the_barrier_before_maturity);
        if(the_process_hits_the_barrier_before_maturity){
            running_sum_for_monte_carlo+=factor*asset_price_process; 
        }
    }
    EXPECT_NEAR(running_sum_for_monte_carlo/num_of_paths,barrier_level,0.99999);
}


TEST(TimeChangedBlackSCholesUpperBarrierMonteCarloTest, CheckTheValueOfApproximation){
    double drift=0.02;
    double volatility=0.2;
    double maturity=1.0;
    double sde_initial_value=100;
    double strike=100; 
    double barrier_level=110;
    unsigned int num_of_subdivisions=10;
    unsigned long int num_of_paths=10000000;

    for(int j=0;j<100;++j){
        double result=euler_maruyama::TimeChangedBlackScholesMonteCarloUpAndInCall(drift,volatility,maturity,sde_initial_value,
                strike,barrier_level,num_of_subdivisions,num_of_paths);
        double expected=black_scholes_explicit_expectation::UpAndInCall(drift,volatility,maturity,strike,sde_initial_value,barrier_level);

        std::cerr<<"[          ] time "<<j+1<<std::endl;
        std::cerr<<"[          ] approximation by time change:"<<result<<std::endl;
        std::cerr<<"[          ] true value:"<<expected<<std::endl;
    }
}

TEST(TimeChangedBlackScholesMonteCarloUpAndInCallVarianceTest, CheckTheValue){
    double drift=0.02;
    double volatility=0.2;
    double maturity=1.0;
    double sde_initial_value=100;
    double strike=100; 
    double barrier_level=110;
    unsigned int num_of_subdivisions=10;
    unsigned long int num_of_paths=10000000;

    double mean_of_variance=0;
    double variance_of_variance=0;

    
    std::cerr<<"[          ] drift:"<<drift<<std::endl;
    std::cerr<<"[          ] volatility:"<<volatility<<std::endl;
    std::cerr<<"[          ] maturity:"<<maturity<<std::endl;
    std::cerr<<"[          ] the initial value of the asset price process:"<<sde_initial_value<<std::endl;
    std::cerr<<"[          ] strike:"<<strike<<std::endl;
    std::cerr<<"[          ] level of the barrier:"<<barrier_level<<std::endl;
    std::cerr<<"[          ] number of subdivisions:"<<num_of_subdivisions<<std::endl;
    std::cerr<<"-------results come up from here-------"<<std::endl;
    for(int j=0;j<10;++j){
        double mean=0;
        double variance=0;
        euler_maruyama::TimeChangedBlackScholesMonteCarloUpAndInCallVariance(drift,volatility,maturity,sde_initial_value,
                strike,barrier_level,num_of_subdivisions,num_of_paths,&mean,&variance);

        std::cerr<<"[          ] time "<<j+1<<std::endl;
        std::cerr<<"[          ] mean:"<<mean<<std::endl;
        std::cerr<<"[          ] variance:"<<variance<<std::endl;

        statistical_functions::MeanVarianceNext(&mean_of_variance,&variance_of_variance,j,variance);

    }
    std::cerr<<"[          ] Summary "<<std::endl;
    std::cerr<<"[          ] averaged variance:"<<mean_of_variance<<std::endl;
}
