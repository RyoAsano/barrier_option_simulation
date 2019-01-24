#include "prob_density_func.h"
#include <cassert>
#include <cmath>
#include "constant.h"
#include "cumulative_dist_func.h"

namespace prob_density_func {

double StdGaussian(double integration_variable){
    return exp(-integration_variable*integration_variable/2.0)/sqrt(2.0*PI);
}

double Gaussian(double mean, double variance, double integration_variable){
    return exp(-(integration_variable-mean)*(integration_variable-mean)/(2.0*variance))/sqrt(2.0*PI*variance);
}

double BrownianMotionFirstHittingTime(double barrier_level, double time){
    assert(time >= 0 && barrier_level >= 0);
    double result = 0;
    if(time!=0){
    result = barrier_level/sqrt(2.0*PI*time*time*time)*exp(-barrier_level*barrier_level/(2.0*time));
    }
    return result;
}

double OneSidedBrownianBridge(double goal_value, double goal_time, double current_time, double integration_variable){
    assert(current_time<goal_time);
    double result=0;
    if(integration_variable<goal_value){
        double density_argument_plus=(integration_variable-goal_value
                +(1.0-current_time/goal_time)*goal_value);
        double density_argument_minus=(integration_variable-goal_value
                -(1.0-current_time/goal_time)*goal_value);
        double variance=(goal_time-current_time)*current_time/goal_time;
 
        result=(goal_value-integration_variable)/goal_value * goal_time/(goal_time-current_time)*
            (prob_density_func::Gaussian(0,variance,density_argument_plus)
             -prob_density_func::Gaussian(0,variance,density_argument_minus));
    }
    return result;
}
}//namespcace PDF


