//
//  gamma_rand_num_generator.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "rand_num_generator.hpp"


double BetaRandNumGenerator(MersenneTwisterUniformRng &unif_gen)
{
    double PI = 3.14159265358979323846;
    double U = unif_gen.next().value;
    
    double arg = PI*U/2.0;
    return std::cos(arg) * std::cos(arg);
}

double ExpRandNumGenerator(MersenneTwisterUniformRng &unif_gen)
{
    double U = unif_gen.next().value;
    return - std::log(U);
}

double GammaRandNumGenerator(MersenneTwisterUniformRng &unif_gen)
{
    double beta = BetaRandNumGenerator(unif_gen);
    double exp = ExpRandNumGenerator(unif_gen);
    
    return beta * exp;
}
