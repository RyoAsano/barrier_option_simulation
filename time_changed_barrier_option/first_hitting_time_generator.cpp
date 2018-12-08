//
//  first_passage_time_generator.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "first_hitting_time_generator.hpp"

double FirstHittingTimeGenerator(MersenneTwisterUniformRng &unif_gen, double barrier_level)
{
    double a = barrier_level;
    double gamma = GammaRandNumGenerator(unif_gen);
    
    return a*a/(2.0*gamma);
}
