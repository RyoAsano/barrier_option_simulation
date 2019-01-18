#pragma once

namespace random_number_generator_func {

double Exponential(double rate, double prob);

double Arcsine(double prob);

double GammaShapeAHalfScaleOne(double unif1, double unif2);

double BrownianMotionFirstHittingTime(double barrier_level, double unif1, double unif2);

}//namespace random_number_generator_func
