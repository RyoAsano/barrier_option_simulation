//
//  function.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/26.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "function.hpp"

using namespace QuantLib;  
using namespace boost::math::quadrature;

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

