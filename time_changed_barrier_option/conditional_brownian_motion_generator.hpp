//
//  conditional_brownian_motion_generator.hpp
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

using namespace QuantLib;
using namespace boost::numeric::ublas;


double lambda(double s, double t, double y, double x);


double IteratedRandomOperators(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                               unsigned long int N, const VectorFields &V, unsigned int j_star,
                               const boost::function<double (vector<double>)> f, double t, double y, vector<double> x);


/*
 conditional_BM_generator generates a Euler--Maruyama approximated path of the conditional Brownian motion given the first hitting time t to the level y.
 N denotes the grid size of the partion for the Euler--Maruyama approximation.
*/
std::vector<double> conditional_BM_generator(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen, unsigned long int N, double t, double y);

/* conditional_brownian_motion_generator_hpp */
