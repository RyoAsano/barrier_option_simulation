#pragma once

namespace test_util {

double GaussianUnusualMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order);

double GaussianTruncatedMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order);

double GaussianSignedMoment(double rv_coeff, double fixed_coeff, double fixed_val, unsigned int moment_order);

double OneSidedBrownianBridgeSignedMoment(double goal_value, double goal_time, double current_time, unsigned int moment_order);
}//namespace test_util
