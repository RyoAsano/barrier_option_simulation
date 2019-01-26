#pragma once

namespace black_scholes_explicit_expectation {

double VanillaCall(double drift, double volatility, double maturity, double strike, double initial_value);

double VanillaPut(double drift, double volatility, double maturity, double strike, double initial_value);

double UpAndInCall(double drift, double volatility, double maturity, double strike, double initial_value, double barrier_level);

double UpAndOutCall(double drift, double volatility, double maturity, double strike, double initial_value, double barrier_level);

}//namespace black_scholes_explicit_expectation
