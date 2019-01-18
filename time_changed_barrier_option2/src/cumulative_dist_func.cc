#include "cumulative_dist_func.h"
#include <cmath>

namespace cumulative_dist_func {

double Gaussian(double prob_less_than){
    return 0.5*erfc(-prob_less_than*sqrt(0.5));
}
}//namespace cumulative_dist_func
