//
//  function.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/26.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "function.hpp"
#include <cassert>

using namespace QuantLib;  
using namespace boost::math::quadrature;
using namespace std;

double ExpInverseCumul(double unif)
{
    return -log(unif);
}

int Lambda1GenFunc(double unif)
{
    if(unif <= 0.5)
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

int Lambda2GenFunc(double unif)
{
    if(unif <= 0.5)
    {
        return 1;
    }
    else
    {
        return -1;
    }
}

double DensityFuncOfBMFirstHittingTime(double barrier_level, double t)
{
    double PI = 3.14159265358979323846;
    double a = barrier_level;

    double result = a/sqrt(2*PI*t*t*t) * exp(-a*a/(2*t));
        return result; 
}



double CDFOfBMFirstHittingTime(double barrier_level, double T)
{
    double PI = 3.14159265358979323846;
    double a = barrier_level;
    double error;
    auto density = [=](double t){return a/sqrt(2*PI*t*t*t) * exp(-a*a/(2*t));};
    
    double result = gauss_kronrod<double, 30>::integrate(density, 0, T, 0, 0, &error);
    
    return result;
}


//cumulative distribution function of the first hitting time of geometric brownian motion (martingale) at an upper barrier.
//P[tau <= T], tau = inf{t>=0 | x exp((mu-sigma^2/2) t + sigma W(t)) >= b}.
double CDFOfGeometricBMFirstHittingTime(double mu, double sigma, double T, double barrier_level, double init_value)
{
    CumulativeNormalDistribution CDF;
    double b = barrier_level;
    double x = init_value;
    
    double h_plus = (-1.0/sigma * log(b/x) + (-mu/sigma + sigma/2.0) * T)/sqrt(T);
    double h_minus = (-1.0/sigma * log(b/x) - (-mu/sigma + sigma/2.0) * T)/sqrt(T);
    
    return CDF(h_minus) + pow(x/b, -2.0*mu/(sigma*sigma)+1.0) * CDF(h_plus);
}


double BlackScholesTrueExpEurCall(double init_value, double drift, double volatility, double maturity, double strike)
{
    double x = init_value;
    double r = drift;
    double sigma = volatility;
    double T = maturity;
    double K = strike;
    CumulativeNormalDistribution CDF;
    
    double delta_plus = 1.0/(sigma*sqrt(T)) * (log(x/K) + (r+sigma*sigma/2.0)*T);
    double delta_minus = 1.0/(sigma*sqrt(T)) * (log(x/K) + (r-sigma*sigma/2.0)*T);
    
    return x*exp(r*T)*CDF(delta_plus) - K*CDF(delta_minus);
}

double BlackScholesTrueExpEurPut(double init_value, double drift, double volatility, double maturity, double strike)
{
    double x = init_value;
    double r = drift;
    double sigma = volatility;
    double T = maturity;
    double K = strike;
    CumulativeNormalDistribution CDF;
    
    double delta_plus = 1.0/(sigma*sqrt(T)) * (log(x/K) + (r+sigma*sigma/2.0)*T);
    double delta_minus = 1.0/(sigma*sqrt(T)) * (log(x/K) + (r-sigma*sigma/2.0)*T);
    
    return K*CDF(-delta_minus) - x *exp(r*T)*CDF(-delta_plus);
}

double BlackScholesTrueExpDigCall(double init_value, double drift, double volatility, double maturity, double strike)
{
    double x = init_value;
    double r = drift;
    double sigma = volatility;
    double T = maturity;
    double K = strike;
    CumulativeNormalDistribution CDF;
    
    double delta_minus = 1.0/(sigma*sqrt(T)) * (log(x/K) + (r-sigma*sigma/2.0)*T);
    
    return CDF(delta_minus);

}

double BlackScholesTrueExpDigPut(double init_value, double drift, double volatility, double maturity, double strike)
{
    double x = init_value;
    double r = drift;
    double sigma = volatility;
    double T = maturity;
    double K = strike;
    CumulativeNormalDistribution CDF;

    double delta_minus = 1.0/(sigma*sqrt(T)) * (log(x/K) + (r-sigma*sigma/2.0)*T);

    return CDF(-delta_minus);
    
}


double BlackScholesTrueExpEurCallUpAndOut(double init_value, double drift, double volatility, double maturity, double strike, double barrier_level)
{
    double S = init_value;
    double r = drift;
    double sigma = volatility;
    double T = maturity;
    double E = strike;
    double B = barrier_level;
    CumulativeNormalDistribution CDF;
    
    double alpha = 0.5 * (1.0 - 2.0*r/(sigma*sigma));

    double result = 0;
    
    result += BlackScholesTrueExpEurCall(S, r, sigma, T, E) - BlackScholesTrueExpEurCall(S, r, sigma, T, B) - (B-E)* BlackScholesTrueExpDigCall(S, r, sigma, T, B);
    
    result += - pow(S/B, 2.0*alpha)
    *(BlackScholesTrueExpEurCall(B*B/S, r, sigma, T, E) - BlackScholesTrueExpEurCall(B*B/S, r, sigma, T, B) - (B-E)*BlackScholesTrueExpDigCall(B*B/S, r, sigma, T, B));
    
    return result;
    
}

double BlackScholesTrueExpEurPutUpAndOut(double init_value, double drift, double volatility, double maturity, double strike, double barrier_level)
{
    double S = init_value;
    double r = drift;
    double sigma = volatility;
    double T = maturity;
    double E = strike;
    double B = barrier_level;
    CumulativeNormalDistribution CDF;
    
    double alpha = 0.5 * (1.0 - 2.0*r/(sigma*sigma));
    
    double result = 0;
    
    result = BlackScholesTrueExpEurPut(S, r, sigma, T, E) - pow(S/B, 2.0*alpha)*BlackScholesTrueExpEurPut(B*B/S, r, sigma, T, E);

    return result;
}

double BlackScholesTrueExpEurCallUpAndIn(double init_value, double drift, double volatility, double maturity, double strike, double barrier_level)
{
    double result = BlackScholesTrueExpEurCall(init_value, drift, volatility, maturity, strike)
    -BlackScholesTrueExpEurCallUpAndOut(init_value, drift, volatility, maturity, strike, barrier_level);
    
    return result;
}



double BlackScholesTrueExpEurPutUpAndIn(double init_value, double drift, double volatility, double maturity, double strike, double barrier_level)
{
    double result = BlackScholesTrueExpEurPut(init_value, drift, volatility, maturity, strike)
    - BlackScholesTrueExpEurPutUpAndOut(init_value, drift, volatility, maturity, strike, barrier_level);
    
    return result;
}


double DensityOfConditionalBM(double t, double y, double s, double x){
    double result=0;
    assert(s<t);
    if(s==0){
        result = (y==x)?10000000:0;
    }
    if(x<y){
        double PI = 3.14159265358979323846;
        auto N = [=](double u, double z){return 1.0/sqrt(2.0*PI*u) * exp(-z*z/(2.0*u));};
        double aaa = (t-s)*s/t;
        double bbb = x-y+(1.0-s/t)*y;
        double ccc =x-y-(1.0-s/t)*y;
        result = (y-x)/(t-s)*t/y*(N((t-s)*s/t, x-y+(1.0-s/t)*y)-N((t-s)*s/t, x-y-(1.0-s/t)*y));
    }
    return result;
}

