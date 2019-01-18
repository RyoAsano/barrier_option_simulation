#include  "true_value.h"
#include "constant.h"
#include <cassert>
#include <cmath>


double CDFofNormDist(double prob_less_than){
    return 0.5*erfc(-prob_less_than*sqrt(0.5));
}

double PDFofBMsFirstHittingTime(double barrier_level, double time){
    assert(time >= 0 && barrier_level >= 0);
    double result = 0;
    if(time!=0){
    result = barrier_level/sqrt(2.0*PI*time*time*time)*exp(-barrier_level*barrier_level/(2.0*time));
    }
    return result;
}




