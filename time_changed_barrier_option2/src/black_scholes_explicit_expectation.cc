#include "black_scholes_explicit_expectation.h"
#include "cumulative_dist_func.h"
#include <cmath>
#include <cassert>

namespace black_scholes_explicit_expectation {

double VanillaCall(double drift, double volatility, double maturity, double strike, double initial_value){
    assert(volatility>0 && maturity > 0);
    double delta_plus=1.0/(volatility*sqrt(maturity))*(log(initial_value/strike)+(drift+volatility*volatility/2.0)*maturity);
    double delta_minus = 1.0/(volatility*sqrt(maturity)) * (log(initial_value/strike)+(drift-volatility*volatility/2.0)*maturity);
    return initial_value*exp(drift*maturity)*cumulative_dist_func::StdGaussian(delta_plus) - strike*cumulative_dist_func::StdGaussian(delta_minus);        
}

double VanillaPut(double drift, double volatility, double maturity, double strike, double initial_value){
    assert(volatility>0 && maturity >0);
    double delta_plus=1.0/(volatility*sqrt(maturity))*(log(initial_value/strike)+(drift+volatility*volatility/2.0)*maturity);
    double delta_minus = 1.0/(volatility*sqrt(maturity)) * (log(initial_value/strike)+(drift-volatility*volatility/2.0)*maturity);
    return strike*cumulative_dist_func::StdGaussian(-delta_minus) - initial_value*exp(drift*maturity)*cumulative_dist_func::StdGaussian(-delta_plus);
}

/*
 * the following equations with the notations are given by John C.Hull's book "OPTIONS FUTURES & OTHER DERIVATIVES"".
 * You can find the formula in 18.1 "Types of Exotic Options", p.463.
 */
double UpAndInCall(double drift, double volatility, double maturity, double strike, double initial_value, double barrier_level){
    assert(initial_value<=barrier_level&&"this is UP-and-in option.");
    assert(maturity>0);
    assert(volatility>0);

    double lambda=(drift+volatility*volatility/2.0)/(volatility*volatility);
    double vol_times_sqrt_maturity=volatility*sqrt(maturity);
    double y=log(barrier_level*barrier_level/(initial_value*strike))/vol_times_sqrt_maturity
                +lambda*vol_times_sqrt_maturity;
    double x1=log(initial_value/barrier_level)/vol_times_sqrt_maturity
                +lambda*vol_times_sqrt_maturity;
    double y1=log(barrier_level/initial_value)/vol_times_sqrt_maturity
                +lambda*vol_times_sqrt_maturity;
    return initial_value*exp(drift*maturity) *cumulative_dist_func::StdGaussian(x1)
            -strike*cumulative_dist_func::StdGaussian(x1-vol_times_sqrt_maturity)
            -initial_value*exp(drift*maturity)*pow(barrier_level/initial_value,2.0*lambda)
                        *(cumulative_dist_func::StdGaussian(-y)-cumulative_dist_func::StdGaussian(-y1))
            +strike*pow(barrier_level/initial_value,2.0*lambda-2.0)
                        *(cumulative_dist_func::StdGaussian(-y+vol_times_sqrt_maturity)-cumulative_dist_func::StdGaussian(-y1+vol_times_sqrt_maturity));
}

double UpAndOutCall(double drift, double volatility, double maturity, double strike, double initial_value, double barrier_level){
    return black_scholes_explicit_expectation::VanillaCall(drift,volatility,maturity,strike,initial_value)
            -black_scholes_explicit_expectation::UpAndInCall(drift,volatility,maturity,strike,initial_value,barrier_level);
}
}//namespace black_sholes_explicit_expectation
