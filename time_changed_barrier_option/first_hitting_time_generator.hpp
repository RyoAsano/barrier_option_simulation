//
//  first_passage_time_generator.hpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//


#pragma once
#include <ql/quantlib.hpp>
#include <ql/methods/montecarlo/sample.hpp>
#include <cmath>
#include "rand_num_generator.hpp"

using namespace QuantLib;

//generates a random number of the first hitting time of a 1-dim brownian motion to the barrier_level.
double FirstHittingTimeGenerator(MersenneTwisterUniformRng &unif_gen, double barrier_level);

 /* first_passage_time_generator_hpp */
