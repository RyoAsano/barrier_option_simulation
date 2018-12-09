//
//  generator_Euler_Maruyama_merger.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/09.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "generator_Euler_Maruyama_merger.hpp"

using namespace boost::numeric::ublas;
using namespace QuantLib;

vector<double> EulerMaruyamaSchemeMerger(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                         const MersenneTwisterUniformRng &unif_gen,
                                         unsigned long int N_EM_cond, unsigned long int N_EM, const VectorFields &V,
                                         unsigned int j_star, double T, double barrier_level, vector<double> x)
{
    double tau = GeneratorFirstHittingTime(unif_gen, barrier_level);
    
    vector<double> running_x = x;
    
    running_x = EulerMaruyamaSchemeWithConditionalBM(norm_rand_gen, N_EM_cond, V, j_star, tau, barrier_level, running_x);
    running_x = EulerMaruyamaScheme(norm_rand_gen, N_EM, V, T, running_x);
    
    return running_x;
}

double IteratedRandomOperatorForEulerMaruyamaMerger(const BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                                    const MersenneTwisterUniformRng &unif_gen,
                                                    unsigned long int N_EM_cond, unsigned long int N_EM, const VectorFields &V, unsigned int j_star, double T, double barrier_level,
                                                    const boost::function<double (vector<double>)> f, vector<double> x)
{
    vector<double> running_x = EulerMaruyamaSchemeMerger(norm_rand_gen, unif_gen, N_EM_cond, N_EM, V, j_star, T, barrier_level, x);
    
    return f(running_x);
}
