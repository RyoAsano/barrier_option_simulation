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

/* function_hpp */
