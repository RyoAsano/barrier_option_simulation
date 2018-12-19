//
//  generator_Euler_Maruyama_merger.hpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/09.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#pragma once

#include <iostream>
#include <limits>
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
#include "generator_Euler_Maruyama.hpp"
#include "generator_Euler_Maruyama_conditional.hpp"
#include "generator_first_hitting_time.hpp"


//simply merges the time-changed and normal SDEs (without stopping condition).
boost::numeric::ublas::vector<double> EulerMaruyamaSchemeMerger(const QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                                                const QuantLib::MersenneTwisterUniformRng &unif_gen,
                                                                unsigned long int N_EM_cond, unsigned long int N_EM, const VectorFields &V_time_changed, const VectorFields &V_normal,
                                                                unsigned int j_star, double T, double barrier_level, boost::numeric::ublas::vector<double> x);

//entails the stopping condition.
double IteratedRandomOperatorForEulerMaruyamaMerger(const QuantLib::BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                                    const QuantLib::MersenneTwisterUniformRng &unif_gen,
                                                    unsigned long int N_EM_cond, unsigned long int N_EM, const VectorFieldsTimeChanged &V_time_changed, const VectorFields &V_normal,
                                                    unsigned int j_star, double T, double barrier_level,
                                                    const boost::function<double (boost::numeric::ublas::vector<double>)> f, boost::numeric::ublas::vector<double> x);

/* generator_Euler_Maruyama_merger_hpp */
