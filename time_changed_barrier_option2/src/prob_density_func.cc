#include "prob_density_func.h"
#include <cassert>
#include <cmath>
#include "constant.h"

namespace prob_density_func {

double StdGaussian(double arg){
    return 1.0/sqrt(2.0*PI)*exp(-arg*arg/2.0);
}

double BrownianMotionFirstHittingTime(double barrier_level, double time){
    assert(time >= 0 && barrier_level >= 0);
    double result = 0;
    if(time!=0){
    result = barrier_level/sqrt(2.0*PI*time*time*time)*exp(-barrier_level*barrier_level/(2.0*time));
    }
    return result;
}
}//namespcace PDF


