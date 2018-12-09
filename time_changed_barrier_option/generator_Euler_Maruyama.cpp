//
//  generator_Euler_Maruyama.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/09.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "generator_Euler_Maruyama.hpp"

using namespace QuantLib;
using namespace boost::numeric::ublas;

vector<double> EulerMaruyamaScheme(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                   unsigned long int N, const VectorFields &V, double T, vector<double> x)
{
    vector<double> running_x = x;           //running_x takes on the first argument of the operator Qf that is supposed to update recursively.
    unsigned int d = V.GetNumOfVecFields();
    
    //impliment the recursive argument up to k <= N-1.
    for(unsigned long int i=0; i<N; i++)
    {
        vector<double> running_x_before_update = running_x;
        running_x = running_x_before_update;

        for(unsigned int j=0; j <=d; j++)
        {
            if(j==0)
            {
                running_x += V.GetVal(0, running_x_before_update) * (T/N);
            }
            else
            {
                double Z = norm_rand_gen.next().value;
                running_x += V.GetVal(j, running_x_before_update) * sqrt(T/N) * Z;
            }
        }
    }
    
    return running_x;
    
}


double IteratedRandomOperatorForEulerMaruyama(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                              unsigned long int N, const VectorFields &V, const boost::function<double (vector<double>)> f, double T, vector<double> x)
{
    vector<double> running_x = EulerMaruyamaScheme(norm_rand_gen, N, V, T, x);

    return f(running_x);
    
}
