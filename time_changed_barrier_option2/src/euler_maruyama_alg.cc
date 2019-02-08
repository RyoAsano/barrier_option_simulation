#include "euler_maruyama_alg.h"
#include <ql/quantlib.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include <cmath>
#include <unistd.h>
#include <time.h>
#include "random_number_generator_func.h"
#include "black_scholes_explicit_expectation.h"
#include "statistical_functions.h"

namespace euler_maruyama {

void BlackScholesOneStep(double drift, double volatility, double std_gaussian_rand_num, double time_step, double* running_process){
    double spot = *running_process;
    *running_process += drift*spot*time_step + volatility*spot*sqrt(time_step)*std_gaussian_rand_num;
}

void TimeChangedBlackScholesUpperBarrierOneStep(double drift, double volatility, double stochastic_integrator_increment, 
        double time_increment, double *asset_price_process_ptr, double *barrier_hitting_time_ptr, double *log_density_process_ptr){
    assert(volatility>0);
    assert((*asset_price_process_ptr)>0);
    *asset_price_process_ptr+=stochastic_integrator_increment;
    *barrier_hitting_time_ptr+=time_increment/(volatility*volatility*(*asset_price_process_ptr)*(*asset_price_process_ptr));
    *log_density_process_ptr+=-drift*drift*time_increment/(2.0*volatility*volatility*volatility*volatility*
                                                                (*asset_price_process_ptr)*(*asset_price_process_ptr))
                              +drift*stochastic_integrator_increment/(volatility*volatility*(*asset_price_process_ptr));
}

void TimeChangedBlackScholesUpperBarrierOnePath(double drift, double volatility, double maturity, double sde_initial_value,
        double barrier_level, double brownian_bridge_goal_time, unsigned int num_of_subdivisions, const double *norm_rand_num_array,
        double *asset_price_process_ptr, double *barrier_hitting_time_ptr, double *log_density_process_ptr, 
        double *factor_ptr, bool *the_process_hits_the_barrier_before_maturity_ptr){
    //assert(num_of_subdivisions==(sizeof(norm_rand_num_array)/sizeof(norm_rand_num_array[0]))&&
     //       "the array of std. norm. rand. nums should be of the num_of_subdivisions number of elements.");
    assert(barrier_level>sde_initial_value);

    //initialization.
    *asset_price_process_ptr=sde_initial_value;
    *factor_ptr=1.0;
    *barrier_hitting_time_ptr=0;
    *log_density_process_ptr=0;
    *the_process_hits_the_barrier_before_maturity_ptr=true;

    double brownian_bridge_goal_value=barrier_level-sde_initial_value; //which is \phi(x)=b-x
    double running_one_sided_brownian_bridge=0;
    double time_increment=brownian_bridge_goal_time/num_of_subdivisions;
    double current_time=0;
    for(int k=0;k<num_of_subdivisions;++k){
        double one_sided_brownian_bridge_increment=(k==num_of_subdivisions-1)?
                (brownian_bridge_goal_value-running_one_sided_brownian_bridge):
                random_number_generator_func::OneSidedBrownianBridgeIncrement(brownian_bridge_goal_value,brownian_bridge_goal_time,
                        current_time,running_one_sided_brownian_bridge,time_increment,norm_rand_num_array[k],factor_ptr);
        euler_maruyama::TimeChangedBlackScholesUpperBarrierOneStep(drift,volatility,one_sided_brownian_bridge_increment,
                time_increment,asset_price_process_ptr,barrier_hitting_time_ptr,log_density_process_ptr);
        if((*asset_price_process_ptr)<=0 || (*barrier_hitting_time_ptr)>maturity){
            *the_process_hits_the_barrier_before_maturity_ptr=false;
            break;
        }
        running_one_sided_brownian_bridge+=one_sided_brownian_bridge_increment;
        current_time+=time_increment;
    }
}

double TimeChangedBlackScholesMonteCarloUpAndInCall(double drift, double volatility, double maturity, double sde_initial_value,
        double strike, double barrier_level, unsigned int num_of_subdivisions, unsigned long num_of_paths){
    //time_t now=time(0);
    //unsigned long seed=getpid()+(unsigned long)now;

    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    typedef QuantLib::MersenneTwisterUniformRng MersenneTwister;
    typedef QuantLib::BoxMullerGaussianRng<MersenneTwister> BoxMuller;
    QuantLib::RandomSequenceGenerator<MersenneTwister> unif_gen(2,MersenneTwister(seed));
    QuantLib::RandomSequenceGenerator<BoxMuller> norm_rand_num_array_gen(num_of_subdivisions,BoxMuller(MersenneTwister(seed)));

    double running_mean_for_monte_carlo=0;
    for(int i=0;i<num_of_paths;++i){
        std::vector<double> unif_vector=unif_gen.nextSequence().value;
        double brownian_bridge_goal_time=random_number_generator_func::BrownianMotionFirstHittingTime(
                barrier_level-sde_initial_value,unif_vector[0],unif_vector[1]);
        double asset_price_process;
        double barrier_hitting_time;
        double log_density_process;
        double factor;
        bool the_process_hits_the_barrier_before_maturity;
        euler_maruyama::TimeChangedBlackScholesUpperBarrierOnePath(drift,volatility,maturity,sde_initial_value,
                barrier_level,brownian_bridge_goal_time,num_of_subdivisions,&(norm_rand_num_array_gen.nextSequence().value[0]),
                &asset_price_process,&barrier_hitting_time,&log_density_process,&factor,&the_process_hits_the_barrier_before_maturity);
        double sample=0;
        if(the_process_hits_the_barrier_before_maturity){
            sample=factor*exp(log_density_process)*
                            black_scholes_explicit_expectation::VanillaCall(drift,volatility,maturity-barrier_hitting_time,strike,barrier_level);
        }
        running_mean_for_monte_carlo=(double)i/(double)(i+1)*running_mean_for_monte_carlo+sample/(double)(i+1);
    }
    return running_mean_for_monte_carlo;
}

void TimeChangedBlackScholesMonteCarloUpAndInCallVariance(double drift, double volatility, double maturity, double sde_initial_value,
        double strike, double barrier_level, unsigned int num_of_subdivisions, unsigned long num_of_paths, double *mean_ptr, double *variance_ptr){

    QuantLib::BigInteger seed = QuantLib::SeedGenerator::instance().get();
    typedef QuantLib::MersenneTwisterUniformRng MersenneTwister;
    typedef QuantLib::BoxMullerGaussianRng<MersenneTwister> BoxMuller;
    QuantLib::RandomSequenceGenerator<MersenneTwister> unif_gen(2,MersenneTwister(seed));
    QuantLib::RandomSequenceGenerator<BoxMuller> norm_rand_num_array_gen(num_of_subdivisions,BoxMuller(MersenneTwister(seed)));

    *mean_ptr=0;
    *variance_ptr=0;
    for(int i=0;i<num_of_paths;++i){
        double next_sample=0;
        std::vector<double> unif_vector=unif_gen.nextSequence().value;
        double brownian_bridge_goal_time=random_number_generator_func::BrownianMotionFirstHittingTime(
                barrier_level-sde_initial_value,unif_vector[0],unif_vector[1]);
        double asset_price_process;
        double barrier_hitting_time;
        double log_density_process;
        double factor;
        bool the_process_hits_the_barrier_before_maturity;
        euler_maruyama::TimeChangedBlackScholesUpperBarrierOnePath(drift,volatility,maturity,sde_initial_value,
                barrier_level,brownian_bridge_goal_time,num_of_subdivisions,&(norm_rand_num_array_gen.nextSequence().value[0]),
                &asset_price_process,&barrier_hitting_time,&log_density_process,&factor,&the_process_hits_the_barrier_before_maturity);
        if(the_process_hits_the_barrier_before_maturity){
            next_sample=factor*exp(log_density_process)*
                black_scholes_explicit_expectation::VanillaCall(drift,volatility,maturity-barrier_hitting_time,strike,barrier_level);
        }
        statistical_functions::MeanVarianceNext(mean_ptr,variance_ptr,i,next_sample);
    }
}
}//namespace euler_maruyama
