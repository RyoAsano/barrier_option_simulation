#pragma once

namespace statistical_functions {

void MeanVarianceNext(double *running_mean_ptr, double *running_variance_ptr, unsigned int current_sample_size, double next_sample_value);

enum AcceptanceLevel {NinetyFive, OneOverAHundred, OneOverAThousand};

unsigned long long SampleSizeGenerator(double accuracy, double variance, statistical_functions::AcceptanceLevel level);

}//namespace statistical_functions
