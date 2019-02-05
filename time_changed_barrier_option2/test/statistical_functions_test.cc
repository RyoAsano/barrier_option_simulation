#include <gtest/gtest.h>
#include "statistical_functions.h"
#include "cumulative_dist_func.h"
#include <cmath>

TEST(MeanVarianceNextCheck, ExactValue){
    double sample_array[]={10,20,30};
    unsigned int running_sample_size=1;
    unsigned int sample_size=sizeof(sample_array)/sizeof(sample_array[0]);
    double running_mean=0;
    double running_variance=0;

   for(int i=0;i<sample_size;++i){
       statistical_functions::MeanVarianceNext(&running_mean,&running_variance,i,sample_array[i]);
   }
   EXPECT_DOUBLE_EQ(running_mean,20);
   EXPECT_DOUBLE_EQ(running_variance,200/3.0);
}

TEST(SampleSize, CheckProbability){
    double variance=1.0;
    double accuracy=0.1;

    unsigned long n_ninty_five=statistical_functions::SampleSizeGenerator(accuracy,variance,statistical_functions::NinetyFive);
    unsigned long n_one_over_a_hundred=statistical_functions::SampleSizeGenerator(accuracy,variance,statistical_functions::OneOverAHundred);
    unsigned long n_one_over_a_thousand=statistical_functions::SampleSizeGenerator(accuracy,variance,statistical_functions::OneOverAThousand);

    std::cerr<<"[        ] variance:"<<variance <<" accuracy:"<<accuracy<<std::endl;
    std::cerr<<"[        ] sample size for 95%:"<<n_ninty_five<<std::endl;
    std::cerr<<"[        ] sample size for 99%:"<<n_one_over_a_hundred<<std::endl;
    std::cerr<<"[        ] sample size for 99.9%:"<<n_one_over_a_thousand<<std::endl;
    EXPECT_NEAR(2.0*cumulative_dist_func::StdGaussian(accuracy*sqrt(n_ninty_five/variance))-1.0,0.95, 0.009);
    EXPECT_NEAR(2.0*cumulative_dist_func::StdGaussian(accuracy*sqrt(n_one_over_a_hundred/variance))-1.0,0.99, 0.009);
    EXPECT_NEAR(2.0*cumulative_dist_func::StdGaussian(accuracy*sqrt(n_one_over_a_thousand/variance))-1.0,0.999, 0.0009);
}
