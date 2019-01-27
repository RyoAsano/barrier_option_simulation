#include "euler_maruyama_alg.h"
#include <ql/quantlib.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include <cmath>
#include <unistd.h>
#include <time.h>
#include "random_number_generator_func.h"

namespace euler_maruyama {

void BlackScholesOneStep(double drift, double volatility, double std_gaussian_rand_num, double time_step, double* running_process){
    double spot = *running_process;
    *running_process += drift*spot*time_step + volatility*spot*sqrt(time_step)*std_gaussian_rand_num;
}

void TimeChangedBlackScholesUpperBarrierOneStep(double volatility, double stochastic_integrator_increment, 
        double time_increment, double *asset_price_process_ptr, double *barrier_hitting_time_ptr){
    assert(volatility>0);
    assert((*asset_price_process_ptr)>0);
    *asset_price_process_ptr+=stochastic_integrator_increment;
    *barrier_hitting_time_ptr+=time_increment/(volatility*volatility*(*asset_price_process_ptr)*(*asset_price_process_ptr));
}

void TimeChangedBlackScholesUpperBarrierOnePath(double volatility, double maturity, double sde_initial_value,
        double barrier_level, double brownian_bridge_goal_time, unsigned int num_of_subdivisions, const double *norm_rand_num_array,
        double *asset_price_process_ptr, double *barrier_hitting_time_ptr, bool *the_process_hits_the_barrier_before_maturity_ptr){
    assert(num_of_subdivisions==(sizeof(norm_rand_num_array)/sizeof(norm_rand_num_array[0]))&&
            "the array of std. norm. rand. nums should be of the num_of_subdivisions number of elements.");
    assert(barrier_level>sde_initial_value);

    //initialization.
    *asset_price_process_ptr=sde_initial_value;
    *barrier_hitting_time_ptr=0;
    *the_process_hits_the_barrier_before_maturity_ptr=true;

    double brownian_bridge_goal_value=barrier_level-sde_initial_value; //which is \phi(x)=b-x
    double running_factor=1.0;
    double running_one_sided_brownian_bridge=0;
    double time_increment=brownian_bridge_goal_time/num_of_subdivisions;
    double current_time=0;
    for(int k=0;k<num_of_subdivisions;++k){
        double one_sided_brownian_bridge_increment=(k==num_of_subdivisions-1)?
                (brownian_bridge_goal_value-running_one_sided_brownian_bridge):
                random_number_generator_func::OneSidedBrownianBridgeIncrement(brownian_bridge_goal_value,brownian_bridge_goal_time,
                        current_time,running_one_sided_brownian_bridge,time_increment,norm_rand_num_array[k],&running_factor);
        euler_maruyama::TimeChangedBlackScholesUpperBarrierOneStep(volatility,one_sided_brownian_bridge_increment,
                time_increment,asset_price_process_ptr,barrier_hitting_time_ptr);
        if((*barrier_hitting_time_ptr)>maturity){
            *the_process_hits_the_barrier_before_maturity_ptr=false;
            break;
        }
        running_one_sided_brownian_bridge+=one_sided_brownian_bridge_increment;
        current_time+=time_increment;
    }
}

void TimeChangedBlackScholesUpperBarrierSample(double volatility, double maturity, double sde_initial_value,
        double borwnian_bridge_goal_value, double brownian_bridge_goal_time, unsigned int num_of_subdivisions, unsigned long num_of_paths){
    time_t now=time(0);
    unsigned long seed=getpid()+(unsigned long)now;
    QuantLib::MersenneTwisterUniformRng unif_gen(seed);
    QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> norm_rand_gen(unif_gen);
}
}//namespace euler_maruyama
