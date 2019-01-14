//
//  generator_Euler_Maruyama.hpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/09.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#pragma once

#include <iostream>
#include <boost/function.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <ql/quantlib.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include "vector_field_class.hpp"
#include <fstream>


boost::numeric::ublas::vector<double> EulerMaruyamaScheme(const QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                                          unsigned long int N, const VectorFields &V, double T, boost::numeric::ublas::vector<double> x);

//returns X(tau \wedge T).
boost::numeric::ublas::vector<double> EulerMaruyamaSchemeWithStoppingCond(const QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                                                          unsigned long int N, const VectorFields &V, double T, boost::numeric::ublas::vector<double> x, double barrier_level, bool *the_process_hits_the_barrier);

//returns f(X(T,x)).
double IteratedRandomOperatorForEulerMaruyama(const QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                              unsigned long int N, const VectorFields &V, const boost::function<double (boost::numeric::ublas::vector<double>)> f, double T, boost::numeric::ublas::vector<double> x);


/* generator_Euler_Maruyama_hpp */
