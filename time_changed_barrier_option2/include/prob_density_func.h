#pragma once

namespace prob_density_func {

double StdGaussian(double integration_variable);

double Gaussian(double mean, double variance, double integration_variable);

double BrownianMotionFirstHittingTime(double barrier_level, double time);

double OneSidedBrownianBridge(double goal_value, double goal_time, double current_time, double integration_variable);
}//namespace pdf
