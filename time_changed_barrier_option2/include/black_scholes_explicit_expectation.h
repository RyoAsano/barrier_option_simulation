#pragma once

namespace black_scholes_explicit_expectation {

double VanillaCall(double drift, double volatility, double maturity, double strike, double initial_value);

double VanillaPut(double drift, double volatility, double maturity, double strike, double initial_value);

}//namespace black_scholes_explicit_expectation
