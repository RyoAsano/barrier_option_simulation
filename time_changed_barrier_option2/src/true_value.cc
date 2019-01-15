#include  "true_value.h"
#include "constant.h"
#include <cmath>

double PDFofFirstHittingTime(double barrier_level, double time){
    return barrier_level/sqrt(2.0*PI*time*time*time)*exp(-barrier_level*barrier_level/(2.0*time));
}


