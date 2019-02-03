#include <gtest/gtest.h>
#include "statistical_functions.h"

TEST(StatisticalFunctionsTest, MeanVarianceCheck){
    double sample_array[]={1.0,2.0,3.0};
    unsigned int running_sample_size=1;
    unsigned int sample_size=sizeof(sample_array)/sizeof(sample_array[0]);
    double running_mean=0;
    double running_variance=0;

   for(int i=0;i<sample_size;++i){
       statistical_functions::mean_variance_next(&running_mean,&running_variance,i,sample_array[i]);
   }
   EXPECT_DOUBLE_EQ(running_mean,2.0);
   EXPECT_DOUBLE_EQ(running_variance,2.0/3.0);
}
