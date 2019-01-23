#pragma once

namespace random_number_generator_func {

double Exponential(double rate, double prob);

double Arcsine(double prob);

double GammaShapeAHalfScaleOne(double unif1, double unif2);

double BrownianMotionFirstHittingTime(double barrier_level, double unif1, double unif2);

void OneSidedBrownianBridgeFromOrigin(double *factor_storage, double *arg_storage, 
        double goal_value, double goal_time, double current_time, double std_norm_rand_num);

double OneSidedBrownianBridgeIncrement(double goal_value, double goal_time,
        double current_time, double current_value, double time_increment, double std_norm_rand_num, double *running_factor_ptr);
}//namespace random_number_generator_func
