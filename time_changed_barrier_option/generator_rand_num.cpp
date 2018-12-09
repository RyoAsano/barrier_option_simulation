//
//  generator_rand_num.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/08.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "generator_rand_num.hpp"

using namespace QuantLib;

double GeneratorBetaRandNum(const MersenneTwisterUniformRng &unif_gen)
{
    double PI = 3.14159265358979323846;
    double U = unif_gen.next().value;
    
    double arg = PI*U/2.0;
    return std::cos(arg) * std::cos(arg);
}

double GeneratorExpRandNum(const MersenneTwisterUniformRng &unif_gen)
{
    double U = unif_gen.next().value;
    return - std::log(U);
}

double GeneratorGammaRandNum(const MersenneTwisterUniformRng &unif_gen)
{
    double beta = GeneratorBetaRandNum(unif_gen);
    double exp = GeneratorExpRandNum(unif_gen);
    
    return beta * exp;
}
