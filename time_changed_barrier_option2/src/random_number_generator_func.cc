#include "random_number_generator_func.h"
#include "constant.h"
#include <cassert>
#include <cmath>

namespace random_number_generator_func {

double Exponential(double rate, double prob){
    assert(rate>0 && 0<=prob && prob<=1);
    return -1.0/rate * log(prob);
}

double Arcsine(double prob){
    assert(0<=prob && prob<=1.0);
    return cos(PI*prob/2.0)*cos(PI*prob/2.0);
}

/*
 * we adopt for generating Gamma(0.5,1.0)-rand. num. the method of
 * Y*Z, where Y and Z follow arcsine and Exp(1.0) dist.'s, resp.
 */
double GammaShapeAHalfScaleOne(double unif1, double unif2){
    assert(0<=unif1 && unif1<=1.0);
    assert(0<=unif2 && unif2<=1.0);
    double arcsine_rand_num=random_number_generator_func::Arcsine(unif1);
    double exp_rand_num=random_number_generator_func::Exponential(1.0, unif2);
    return arcsine_rand_num*exp_rand_num;
}

double BrownianMotionFirstHittingTime(double barrier_level, double unif1, double unif2){
    assert(barrier_level>=0);
    return barrier_level*barrier_level/(2.0*random_number_generator_func::GammaShapeAHalfScaleOne(unif1,unif2));
}


/*
 * this function stores values \Lambda and \lambda in factor_storage and arg_storage, resp.
 * so that we can compute
 * \mathbb{E}[f(\beta_{t,y}(s))]=\mathbb{E}[\Lambda * f(\lambda)],
 * where y, t and s are goal_value, goal_time and current_time, resp.
 */
void OneSidedBrownianBridgeFromOrigin(double *factor_storage, double *arg_storage, 
        double goal_value, double goal_time, double current_time, double std_norm_rand_num){

    assert(goal_value>=0);
    assert(goal_time>=0);  
    assert(current_time>=0);
    assert(goal_time>current_time&&
            "If you want to deal with the case where goal_time==current_time,\
            you should not use the function \"OneSidedBrownianBridgeFromOrigin\".\
            Please separate the case from others for implementation.");

    *factor_storage=1.0-std_norm_rand_num/goal_value*sqrt(current_time*goal_time/(goal_time-current_time));
   
    *arg_storage=goal_value 
                     - abs(std_norm_rand_num*
                                sqrt((goal_time-current_time)*current_time/goal_time)
                           -goal_value*(goal_time-current_time)/goal_time);
}

double OneSidedBrownianBridgeIncrement(double goal_value, double goal_time,
        double current_time, double current_value, double time_increment, double std_norm_rand_num, double *running_factor_ptr){

    double next_time=current_time+time_increment;
    
    assert(goal_value>current_value);
    assert(time_increment>0);
    assert(goal_time>next_time&&
            "If you want to deal with the case where goal_time==next_time,\
            you should not use the function \"OneSidedBrownianBridgeIncrement\".\
            Please separate the case from others for implementation.");
    assert(current_time>=0);

    *running_factor_ptr *=1.0-std_norm_rand_num/(goal_value-current_value)*
                                     sqrt(goal_time-current_time)/(goal_time-next_time)*time_increment;
   
    return goal_value-current_value-abs(
            std_norm_rand_num*sqrt((goal_time-next_time)/(goal_time-current_time)*time_increment)
            -(goal_value-current_value)*(goal_time-next_time)/(goal_time-current_time)
            );
}
}//namespace random_number_generator_func


