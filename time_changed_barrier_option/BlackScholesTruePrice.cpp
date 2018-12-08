#include "BlackScholesTruePrice.h"
#include <cmath>
#include "GaussianDistributionClass.h"


double BlackScholesTruePrice(double InitialValue, double drift, double volatility, double maturity, double strike, double BarrierLevel)
{
    double x = InitialValue;
    double mu = drift;
    double sigma = volatility;
    double T = maturity;
    double K = strike;
    double b = BarrierLevel;
    StandardGaussianDistribution N;

    double result = 0;


    result += x*exp(mu*T)*N.cdf(DeltaPlus(mu, sigma, T, x/b));
    result += -K*N.cdf(DeltaMinus(mu, sigma, T, x/b));
    result += b*pow(x/b, -2*mu/(sigma*sigma))*exp(mu*T)*(N.cdf(DeltaPlus(mu, sigma, T, b*b/(x*K))) - N.cdf(DeltaPlus(mu, sigma, T, b/x)));
    result += -K*pow(x/b, -2*mu/(sigma*sigma)+1)*(N.cdf(DeltaMinus(mu, sigma, T, b*b/(x*K))) - N.cdf(DeltaMinus(mu, sigma, T, b/x)));

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