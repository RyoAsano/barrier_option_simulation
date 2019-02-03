#pragma once

namespace statistical_functions {

void mean_variance_next(double *running_mean_ptr, double *running_variance_ptr, unsigned int current_sample_size, double next_sample_value);

}//namespace statistical_functions
