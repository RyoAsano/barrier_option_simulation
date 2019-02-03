#include "statistical_functions.h"

namespace statistical_functions {
void mean_variance_next(double *running_mean_ptr, double *running_variance_ptr, unsigned int current_sample_size, double next_sample_value){ 
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
}//namespace statistical_functions
