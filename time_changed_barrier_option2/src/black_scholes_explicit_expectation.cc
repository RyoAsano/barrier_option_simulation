#include "black_scholes_explicit_expectation.h"
#include "cumulative_dist_func.h"
#include <cmath>
#include <cassert>

namespace black_scholes_explicit_expectation {

double VanillaCall(double drift, double volatility, double maturity, double strike, double initial_value){
    assert(volatility>0 && maturity > 0);
    double delta_plus=1.0/(volatility*sqrt(maturity))*(log(initial_value/strike)+(drift+volatility*volatility/2.0)*maturity);
    double delta_minus = 1.0/(volatility*sqrt(maturity)) * (log(initial_value/strike)+(drift-volatility*volatility/2.0)*maturity);
    return initial_value*exp(drift*maturity)*cumulative_dist_func::Gaussian(delta_plus) - strike*cumulative_dist_func::Gaussian(delta_minus);        
}

double VanillaPut(double drift, double volatility, double maturity, double strike, double initial_value){
    assert(volatility>0 && maturity >0);
    double delta_plus=1.0/(volatility*sqrt(maturity))*(log(initial_value/strike)+(drift+volatility*volatility/2.0)*maturity);
    double delta_minus = 1.0/(volatility*sqrt(maturity)) * (log(initial_value/strike)+(drift-volatility*volatility/2.0)*maturity);
    return strike*cumulative_dist_func::Gaussian(-delta_minus) - initial_value*exp(drift*maturity)*cumulative_dist_func::Gaussian(-delta_plus);
}
}//namespace black_sholes_explicit_expectation
