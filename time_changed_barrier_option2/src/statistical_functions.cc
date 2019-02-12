#include "statistical_functions.h"
#include <cassert> 

namespace statistical_functions {
void MeanVarianceNext(double *running_mean_ptr, double *running_variance_ptr, unsigned int current_sample_size, double next_sample_value){ 
    if(current_sample_size==0){
        *running_mean_ptr=next_sample_value;
        *running_variance_ptr=0;
    }else{
        double mean_temp=*running_mean_ptr; 
        *running_mean_ptr=(double)current_sample_size/(current_sample_size+1)*(
                *running_mean_ptr+next_sample_value/current_sample_size);
        *running_variance_ptr=(double)current_sample_size/(current_sample_size+1)*(
                *running_variance_ptr+(next_sample_value-mean_temp)*(next_sample_value-mean_temp)/(current_sample_size+1));
    }
}

unsigned long long SampleSizeGenerator(double accuracy, double variance, statistical_functions::AcceptanceLevel level){
    double inverse_of_gaussian_cdf;
    switch(level){
        case AcceptanceLevel::NinetyFive: {
                                              inverse_of_gaussian_cdf=1.95996398454;
                                              break;
                                          }
        case AcceptanceLevel::OneOverAHundred: {
                                                   inverse_of_gaussian_cdf=2.57582930355;
                                                   break;
                                               }
        case AcceptanceLevel::OneOverAThousand: {
                                                    inverse_of_gaussian_cdf=3.29052673149;
                                                    break;
                                                }
        default: {
                    assert(false&&"there is no such an acceptance level defined.");
                 }
    }
    return (unsigned long long)(variance*inverse_of_gaussian_cdf*inverse_of_gaussian_cdf/(accuracy*accuracy));
}
}//namespace statistical_functions
