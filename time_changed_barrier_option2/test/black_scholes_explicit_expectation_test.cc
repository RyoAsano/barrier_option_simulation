#include <gtest/gtest.h>
#include "black_scholes_explicit_expectation.h"
#include <cmath>

TEST(BlackScholesExplicitExpectationTest, ExactValue){
    double drift=0.1;
    double vol=0.2;
    double maturity=0.5;
    double strike=40;
    double init_val=42;

    EXPECT_NEAR(exp(-0.05)*black_scholes_explicit_expectation::VanillaCall(
                drift,vol,maturity,strike,init_val), 4.76, 0.00999999999);

    EXPECT_NEAR(exp(-0.05)*black_scholes_explicit_expectation::VanillaPut(
                drift,vol,maturity,strike,init_val), 0.81, 0.00999999999);
}

