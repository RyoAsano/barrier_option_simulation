#pragma once

namespace test_util {

double GaussianUnusualMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order);

double GaussianTruncatedMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order);

double GaussianSignedMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order);

double OneSidedBrownianBridgeMoment(double goal_value, double goal_time, double current_time, unsigned int moment_order);

void BlackScholesTimeChangedHittingTimeLinearCase(double volatility, double brownian_bridge_initial_value,
        double brownian_bridge_goal_value, double brownian_bridge_goal_time, double current_time,
        double *asset_price_process_ptr, double *barrier_hitting_time_ptr);
}//namespace test_util
