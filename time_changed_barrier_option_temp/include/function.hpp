//
//  function.hpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/26.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#pragma once
#include <ql/quantlib.hpp>
#include <boost/function.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <ql/math/distributions/normaldistribution.hpp>
#include <cmath>


double ExpInverseCumul(double unif);

int Lambda1GenFunc(double unif);

int Lambda2GenFunc(double unif);

double DensityFuncOfBMFirstHittingTime(double barrier_level, double t);

double CDFOfBMFirstHittingTime(double barrier_level, double T);

double CDFOfGeometricBMFirstHittingTime(double mu, double sigma, double T, double barrier_level, double init_value);


//returns the price times exp(rT).
double BlackScholesTrueExpEurCall(double init_value, double drift, double volatility, double maturity, double strike);

//returns the price times exp(rT).
double BlackScholesTrueExpEurPut(double init_value, double drift, double volatility, double maturity, double strike);

double BlackScholesTrueExpDigCall(double init_value, double drift, double volatility, double maturity, double strike);

double BlackScholesTrueExpDigPut(double init_value, double drift, double volatility, double maturity, double strike);

double BlackScholesTrueExpEurCallUpAndOut(double init_value, double drift, double volatility, double maturity, double strike, double barrier_level);

double BlackScholesTrueExpEurPutUpAndOut(double init_value, double drift, double volatility, double maturity, double strike, double barrier_level);

double BlackScholesTrueExpEurCallUpAndIn(double init_value, double drift, double volatility, double maturity, double strike, double barrier_level);

double BlackScholesTrueExpEurPutUpAndIn(double init_value, double drift, double volatility, double maturity, double strike, double barrier_level);

double DensityOfConditionalBM(double t, double y, double s, double x);


/* function_hpp */
