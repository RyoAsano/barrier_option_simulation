//
//  generator_first_hitting_time.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "generator_first_hitting_time.hpp"

double GeneratorFirstHittingTime(MersenneTwisterUniformRng &unif_gen, double barrier_level)
{
    double a = barrier_level;
    double gamma = GeneratorGammaRandNum(unif_gen);
    
    return a*a/(2.0*gamma);
}
