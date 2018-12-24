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
                                         unsigned long int N_EM_cond, unsigned long int N_EM, const VectorFields &V_time_changed, const VectorFields &V_normal,
                                         unsigned int j_star, double T, double barrier_level, vector<double> x)
{
    double tau = GeneratorFirstHittingTime(unif_gen, barrier_level);
    
    vector<double> running_x = x;
    
    running_x = EulerMaruyamaSchemeWithConditionalBM(norm_rand_gen, N_EM_cond, V_time_changed, j_star, tau, barrier_level, running_x);
    running_x.resize(V_normal.GetDimOfStateSpace());
    running_x = EulerMaruyamaScheme(norm_rand_gen, N_EM, V_normal, T, running_x);
    
    return running_x;
}

double IteratedRandomOperatorForEulerMaruyamaMerger(const BoxMullerGaussianRng<QuantLib::MersenneTwisterUniformRng> &norm_rand_gen,
                                                    const MersenneTwisterUniformRng &unif_gen,
                                                    unsigned long int N_EM_cond, unsigned long int N_EM, const VectorFieldsTimeChanged &V_time_changed, const VectorFields &V_normal,
                                                    unsigned int j_star, double T, const boost::function<double (vector<double>)> f, vector<double> x)
{
    double tau_hat = GeneratorFirstHittingTime(unif_gen, abs(V_time_changed.BarrierFunction(x)));
    
    vector<double> running_x = x;
    bool the_process_hits_the_barrier;
    
    running_x = EulerMaruyamaSchemeWithConditionalBMAndStoppingCond(norm_rand_gen, N_EM_cond, V_time_changed, j_star, tau_hat, T, running_x, &the_process_hits_the_barrier);
    
    if(the_process_hits_the_barrier)
    {
        double tau = running_x(V_time_changed.GetBarrierMonitoringIndex());
        running_x.resize(V_normal.GetDimOfStateSpace());
        running_x = EulerMaruyamaScheme(norm_rand_gen, N_EM, V_normal, T-tau, running_x);
        
        return f(running_x);

    }
    else
    {
        return 0;
    }
}
