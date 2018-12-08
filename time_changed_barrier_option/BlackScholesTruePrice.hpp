#pragma once

#include <cmath>
#include <ql/math/distributions/normaldistribution.hpp>

double BlackScholesTruePrice(double InitialValue, double drift, double volatility, double maturity, double strike, double BarrierLevel);
double DeltaPlus(double drift, double volatility, double t, double s);
double DeltaMinus(double drift, double volatility, double t, double s);
