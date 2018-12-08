//
//  operators_for_conditioned_Brownian_motion.hpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/11/21.
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
                               const InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<double(double)>> &exp_rand_gen,
                               const InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<int(double)>> &Lambda1_rand_gen,
                               const InverseCumulativeRsg<RandomSequenceGenerator<MersenneTwisterUniformRng>, boost::function<int(double)>> &Lambda2_rand_gen,
                               unsigned long int N, const VectorFields &V, unsigned int j_star,
                               const boost::function<double (vector<double>)> f, double t, double y, vector<double> x);
/* operators_for_conditioned_Brownian_motion_hpp */

