#include "BlackScholesTruePrice.hpp"

using namespace QuantLib;

double BlackScholesTruePrice(double InitialValue, double drift, double volatility, double maturity, double strike, double BarrierLevel)
{
    double x = InitialValue;
    double mu = drift;
    double sigma = volatility;
    double T = maturity;
    double K = strike;
    double b = BarrierLevel;
    CumulativeNormalDistribution CDF;
    double result = 0;


    result += x*exp(mu*T)*CDF(DeltaPlus(mu, sigma, T, x/b));
    result += -K*CDF(DeltaMinus(mu, sigma, T, x/b));
    result += b*pow(x/b, -2*mu/(sigma*sigma))*exp(mu*T)*(CDF(DeltaPlus(mu, sigma, T, b*b/(x*K))) - CDF(DeltaPlus(mu, sigma, T, b/x)));
    result += -K*pow(x/b, -2*mu/(sigma*sigma)+1)*(CDF(DeltaMinus(mu, sigma, T, b*b/(x*K))) - CDF(DeltaMinus(mu, sigma, T, b/x)));

    return result;
}


double DeltaPlus(double drift, double volatility, double t, double s)
{
    double mu = drift;
    double sigma = volatility;

    return 1/(sigma*sqrt(t)) * (log(s)+(mu+sigma*sigma/2)*t);
}

double DeltaMinus(double drift, double volatility, double t, double s)
{
    double mu = drift;
    double sigma = volatility;

    return 1/(sigma*sqrt(t)) * (log(s)+(mu-sigma*sigma/2)*t);
}
