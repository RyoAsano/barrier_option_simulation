#include <gtest/gtest.h>
#include "true_value.h"

TEST(CDFofNormDist, ProbHalf){
    EXPECT_DOUBLE_EQ(0.5, CDFofNormDist(0));
    EXPECT_LT(CDFofNormDist(-0.05), 0.5);
    EXPECT_LT(0.5, CDFofNormDist(0.05));
}

TEST(CDFofNormDist, Symmetry){
    EXPECT_DOUBLE_EQ(CDFofNormDist(0.5), 1.0-CDFofNormDist(-0.5));
    EXPECT_DOUBLE_EQ(CDFofNormDist(1.0), 1.0-CDFofNormDist(-1.0));
    EXPECT_DOUBLE_EQ(CDFofNormDist(-0.2), 1.0-CDFofNormDist(0.2));
}

TEST(CDFofNormDist, ExactValue){
    EXPECT_NEAR(CDFofNormDist(-2.82), 0.0024, 0.0000999999999); 
    EXPECT_NEAR(CDFofNormDist(-1.22), 0.1112, 0.0000999999999);
    EXPECT_NEAR(CDFofNormDist(1.58), 0.9429, 0.0000999999999);
    EXPECT_NEAR(CDFofNormDist(0.45), 0.6736, 0.0000999999999);
    
}

TEST(PDFofBMsFirstHittingTimeTest, TimeZero){
    EXPECT_DOUBLE_EQ(0,PDFofBMsFirstHittingTime(1, 0));
    EXPECT_DOUBLE_EQ(0,PDFofBMsFirstHittingTime(0,0));
}

TEST(PDFofBMsFirstHittingTimeTest, ExactValue){
    EXPECT_NEAR(PDFofBMsFirstHittingTime(3.0,1.5), 0.03243478222, 0.00000000000999999999);
    EXPECT_NEAR(PDFofBMsFirstHittingTime(2.2,9.5), 0.02323360138, 0.00000000000999999999);
}
