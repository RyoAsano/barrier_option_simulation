//
//  generator_BM.hpp
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

using namespace QuantLib;
using namespace boost::numeric::ublas;


//this is the normal Euler--Maruyama scheme.
double IteratedRandomOperatorForBM(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                   unsigned long int N, const VectorFields &V, const boost::function<double (vector<double>)> f, double T, vector<double> x);


/* generator_BM_hpp */
