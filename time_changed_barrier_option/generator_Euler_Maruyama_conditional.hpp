//
//  generator_Euler_Maruyama_conditional.hpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
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


double lambda(double s, double t, double y, double x);


//returns the realization of the r.v. at time t following the Euler--Maruyama scheme with vector field V.
boost::numeric::ublas::vector<double> EulerMaruyamaSchemeWithConditionalBM(const QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                                                           unsigned long int N, const VectorFields &V, unsigned int j_star, double t, double y, boost::numeric::ublas::vector<double> x);

//returns f(X(t,x)).
double IteratedRandomOperatorForEulerMaruyamaConditional(const QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                                         unsigned long int N, const VectorFields &V, unsigned int j_star,
                                                         const boost::function<double (boost::numeric::ublas::vector<double>)> f, double t, double y, boost::numeric::ublas::vector<double> x);

/*
 GeneratorConditionalBM generates a Euler--Maruyama approximated path of the conditional Brownian motion given the first hitting time t to the level y.
 N denotes the grid size of the partion for the Euler--Maruyama approximation.
*/
std::vector<double> GeneratorConditionalBM(const QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen, unsigned long int N, double t, double y);


/* generator_Euler_Maruyama_conditional_hpp */
