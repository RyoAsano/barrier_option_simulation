//
//  generator_rand_num.hpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#pragma once
#include <ql/quantlib.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include <cmath>


//Generates the Beta(1/2, 1/2) (alpha = 1/2, beta = 1/2) random number from a uniform random variable.
double GeneratorBetaRandNum(QuantLib::MersenneTwisterUniformRng &unif_gen);

//Generates the Exp(1) (lambda = 1) random number.
double GeneratorExpRandNum(QuantLib::MersenneTwisterUniformRng &unif_gen);

//Generates the Gamma(1/2, 1) (alpha = 1/2, lambda = 1) random number.
double GeneratorGammaRandNum(QuantLib::MersenneTwisterUniformRng &unif_gen);

/* generator_rand_num_hpp */
