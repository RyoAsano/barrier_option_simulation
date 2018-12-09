//
//  generator_BM.cpp
//  time_changed_barrier_option
//
//  Created by Asano Ryo on 2018/12/09.
//  Copyright © 2018年 Asano Ryo. All rights reserved.
//

#include "generator_BM.hpp"

using namespace QuantLib;

double IteratedRandomOperatorForBM(const BoxMullerGaussianRng<MersenneTwisterUniformRng> &norm_rand_gen,
                                   unsigned long int N, const VectorFields &V, const boost::function<double (vector<double>)> f, double T, vector<double> x)
{
    vector<double> running_x = x;           //running_x takes on the first argument of the operator Qf that is supposed to update recursively.
    unsigned int d = V.GetNumOfVecFields();
    
    //impliment the recursive argument up to k <= N-1.
    for(unsigned long int i=0; i<N; i++)
    {
        vector<double> running_x_before_update = running_x;
        running_x = running_x_before_update;
        double a = running_x(0);
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

    return f(running_x);
    
}
