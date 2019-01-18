#include <iostream>
#include "gtest/gtest.h"
//#include "true_value.h"
//TEST(CDFofNormDist, ProbHalf){
//    EXPECT_DOUBLE_EQ(0.5, CDFofNormDist(0));
//    EXPECT_LT(0.5, CDFofNormDist(-0.05));
//    EXPECT_LT(CDFofNormDist(0.05), 0);
//}
//
//TEST(CDFofNormDist, Symmetry){
//    EXPECT_DOUBLE_EQ(CDFofNormDist(0.5), 1.0-CDFofNormDist(-0.5));
//    EXPECT_DOUBLE_EQ(CDFofNormDist(1.0), 1.0-CDFofNormDist(-1.0));
//    EXPECT_DOUBLE_EQ(CDFofNormDist(-0.2), 1.0-CDFofNormDist(0.2));
//}
//
//
//TEST(PDFofBMsFirstHittingTimeTest, TimeZero){
//    EXPECT_DOUBLE_EQ(0,PDFofBMsFirstHittingTime(1, 0));
//    EXPECT_DOUBLE_EQ(0,PDFofBMsFirstHittingTime(0,0));
//    EXPECT_NE(0,PDFofBMsFirstHittingTime(0, 0.000000001));
//}
int main(int argc, char **argv){

    ::testing::InitGoogleTest(&argc, argv); 
    return RUN_ALL_TESTS();
}
