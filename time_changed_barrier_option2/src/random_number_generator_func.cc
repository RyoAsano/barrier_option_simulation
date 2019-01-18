#include "random_number_generator_func.h"
#include "constant.h"
#include <cassert>
#include <cmath>

namespace random_number_generator_func {

double Exponential(double rate, double prob){
    assert(rate>0 && prob>=0 && prob<=1);
    return -1.0/rate * log(prob);
}

double Arcsine(double prob){
    assert(prob>=0 && prob<=1.0);
    return cos(PI*prob/2.0)*cos(PI*prob/2.0);
}

/*
 * we adopt for generating Gamma(0.5,1.0)-rand. num. the method of
 * Y*Z, where Y and Z follow arcsine and Exp(1.0) dist.'s, resp.
 */
double GammaShapeAHalfScaleOne(double unif1, double unif2){
    assert(unif1>=0 && unif1<=1.0 && unif2>=0 && unif2<=1.0);
    double arcsine_rand_num=random_number_generator_func::Arcsine(unif1);
    double exp_rand_num=random_number_generator_func::Exponential(1.0, unif2);
    return arcsine_rand_num*exp_rand_num;
}

double BrownianMotionFirstHittingTime(double barrier_level, double unif1, double unif2){
    assert(barrier_level>=0);
    return barrier_level*barrier_level/(2.0*random_number_generator_func::GammaShapeAHalfScaleOne(unif1,unif2));
}
}//namespace random_number_generator_func


